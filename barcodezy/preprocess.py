import sys, os
from skbio import io as skIO
from skbio import DNA as skDNA

def rename_reads(raw_read_file, exp_name) -> None:
    long_reads = [skDNA(seq) for seq in skIO.read(raw_read_file, format='fastq', variant='sanger')]

    for i,longread in enumerate(long_reads):
        # rename the read something more sensible, but keep the old ones too
        longread.metadata['ogid'] = longread.metadata['id']
        longread.metadata['id'] = f'{exp_name}_{i}'

        sys.stdout.write(f'{i+1}/{len(long_reads)}')
        sys.stdout.flush()
        sys.stdout.write('\r')

    sys.stdout.flush()
    sys.stdout.write('\r')

    def data_gen():
        for read in long_reads:
            yield read
    skIO.write(data_gen(), format='fastq', variant='sanger', into=raw_read_file)

def mm2_align(out_dir, exp_name, plasmid_path, trimmed_reads_path) -> None:
    output_bam_file = f'{out_dir}/{exp_name}.bam'
    output_fa_aligned = f'{out_dir}/{exp_name}.aligned.fa'
    output_fa_unaligned = f'{out_dir}/{exp_name}.unaligned.fa.gz'
    rev_complement_cmd = 'while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done'
    if not os.path.exists(output_fa_unaligned):
        # align reads with mm2, convert to bam output
        print(f'Aligning {exp_name} to {os.path.basename(plasmid_path)} using minimap2')
        os.system(f'minimap2 -ax map-ont {plasmid_path} {trimmed_reads_path} \
        --secondary=no -t 3 | samtools view -b -h -@ 2 >{output_bam_file}')
        print(f'Alignment finished. Flipping reverse reads and writing fasta outputs')
        # write positive strand alignments to file
        os.system(f'samtools view -F 2048 -F 4 -F 16 -h {output_bam_file} | samtools fasta - > {output_fa_aligned}')
        # reverse complement negative strand alignments and append to same file
        os.system(f'samtools view -F 2048 -f 16 -h {output_bam_file} | samtools fasta - | {rev_complement_cmd} >> {output_fa_aligned}')
        os.system(f'gzip {output_fa_aligned}') # compress into .fa.gz
        # write unaligned reads to file
        os.system(f'samtools view -f 4 -h {output_bam_file} | samtools fasta - | gzip > {output_fa_unaligned}')