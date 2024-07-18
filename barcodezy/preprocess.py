import sys, os
from skbio import io as skIO
from skbio import DNA as skDNA
from importlib.resources import files
from typing import Dict

def rename_reads(raw_read_file, exp_name) -> None:
    long_reads = [skDNA(seq) for seq in skIO.read(raw_read_file, format='fastq', variant='sanger')]

    for i,longread in enumerate(long_reads):
        # rename the read something more sensible, but keep the old ones too
        longread.metadata['ogid'] = longread.metadata['id']
        longread.metadata['id'] = f'{exp_name}_{i}'

        #sys.stdout.write(f'{i+1}/{len(long_reads)}')
        #sys.stdout.flush()
        #sys.stdout.write('\r')

    #sys.stdout.flush()
    #sys.stdout.write('\r')

    def data_gen():
        for read in long_reads:
            yield read
    print(f'Writing file...')
    skIO.write(data_gen(), format='fastq', variant='sanger', into=raw_read_file)

def mm2_align(out_dir, exp_name, plasmid_path, trimmed_reads_path) -> Dict[str, int]:
    output_bam_file = f'{out_dir}/{exp_name}.bam'
    output_sorted_bam_file = f'{out_dir}/{exp_name}.sorted.bam'
    output_fa_aligned = f'{out_dir}/{exp_name}.aligned.fa'
    output_fa_unaligned = f'{out_dir}/{exp_name}.unaligned.fa.gz'
    rev_complement_cmd = 'while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done'
    # E. coli ref:
    ecoli_ref_fa = files('barcodezy.plasmids').joinpath('Escherichia_coli_gca_001606525.ASM160652v1_.dna_rm.toplevel.fa')
    # Write tmp combined reference of plasmid + ecoli (for quantifying ecoli reads)
    # (without using combined reference, many reads will "align" to both plasmid and ecoli)
    os.system(f'cat {ecoli_ref_fa} {plasmid_path} >{out_dir}/.tmp.ref.fa')
    # Extract plasmid reference name from fasta, for pulling out plasmid reads from the combined alignment
    with open(plasmid_path, 'r') as fh:
        ref_name = str(fh.readline().strip().split(' ')[0][1:])
    if not os.path.exists(output_fa_unaligned):
        # align reads with mm2, convert to bam output
        print(f'\nAligning {exp_name} to {os.path.basename(plasmid_path)} using minimap2')
        os.system(f'minimap2 -ax map-ont {out_dir}/.tmp.ref.fa {trimmed_reads_path} \
        --secondary=no -t 3 | samtools view -b -h -@ 2 >{output_bam_file}')
        # sort and index bam
        os.system(f'samtools sort {output_bam_file} >{output_sorted_bam_file}')
        os.system(f'samtools index {output_sorted_bam_file} && rm {output_bam_file} && rm {out_dir}/.tmp.ref.fa')
        print(f'\nAlignment finished. Flipping reverse reads and writing fasta outputs')

        # write positive strand alignments to fasta (only reads aligned to reference name, not ecoli)
        # SAM flags used:
        ## 2048: supplementary alignments
        ## 4:    unmapped
        ## 16:   reverse strand
        os.system(f'samtools view -F 2048 -F 4 -F 16 -h {output_sorted_bam_file} \'{ref_name}\' | \
                  samtools fasta - > {output_fa_aligned}')
        # reverse complement negative strand alignments and append to same file
        os.system(f'samtools view -F 2048 -f 16 -h {output_sorted_bam_file} \'{ref_name}\' | samtools fasta - | \
                  {rev_complement_cmd} >> {output_fa_aligned}')
        os.system(f'gzip {output_fa_aligned}') # compress into .fa.gz
        # write unaligned reads to file
        os.system(f'samtools view -f 4 -h {output_sorted_bam_file} | samtools fasta - | gzip > {output_fa_unaligned}')
        
        # store counts for aligned reads, unaligned reads, ecoli reads
        summary_dict = {}
        summary_dict['count_ecoli'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                                                    grep -v \'{ref_name}\' | \
                                                    wc -l').read())
        summary_dict['count_aligned'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                                                    grep \'{ref_name}\' | \
                                                    wc -l').read())
        summary_dict['count_unaligned'] = int(os.popen(f'samtools view -f 4 {output_sorted_bam_file} | wc -l').read())

        print(f'\n{summary_dict["count_aligned"]} plasmid reads aligned.')
        print(f'{summary_dict["count_ecoli"]} bacterial reads aligned.')
        print(f'{summary_dict["count_unaligned"]} undetermined reads.\n')
        return summary_dict