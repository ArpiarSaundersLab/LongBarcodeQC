#!/usr/bin/env python
import os
from barcodezy import parser
from barcodezy import preprocess
from barcodezy import barcode_aligner

def main(args=None):
    args = parser.getArgs() # grab user arguments
    
    # inputs
    read_dir = os.path.normpath(args.input)
    ref_plasmid = os.path.normpath(args.plasmid)
    barcode_design = os.path.normpath(args.barcodes)
    insert_length = args.insert_length
    # output
    output_dir = os.path.normpath(args.output)
    experiment_name = os.path.basename(output_dir)
    output_file = f'{output_dir}/{experiment_name}.csv'

    # validate user options 
    parser.validateArgs(read_dir, output_dir, ref_plasmid, barcode_design)

    # remove any adaptor sequences (requires porechop), condense reads into a single file
    trimmed_read_file = f'{output_dir}/{experiment_name}_trimmed.fastq'
    if not os.path.exists(trimmed_read_file):
        cmd_str = (f'porechop -i {read_dir} -o {trimmed_read_file} ' 
                   '--no_split --min_trim_size 25 '
                   '--adapter_threshold 97 --check_reads 2000')
        os.system(cmd_str)
        preprocess.rename_reads(trimmed_read_file, experiment_name)

    # map trimmed reads to plasmid with minimap2
    summary_align_counts = preprocess.mm2_align(output_dir, 
                                                experiment_name, 
                                                ref_plasmid, 
                                                trimmed_read_file)
    # expected path to mm2 aligned reads:
    mapped_reads = f'{output_dir}/{experiment_name}.aligned.fa.gz'

    # generate barcode alignment & read stats written to csv output
    barcode_aligner.barcode_scores(mapped_reads, barcode_design, insert_length, output_file)

if __name__ == "__main__":
    main()