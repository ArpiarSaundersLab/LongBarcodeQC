#!/usr/bin/env python
import os
from barcodezy import parser
from barcodezy import preprocess
from barcodezy import barcode_aligner
from barcodezy import analysis
from importlib.resources import files

def main(args=None):
    args = parser.getArgs() # grab user arguments
    
    # inputs
    read_dir = os.path.normpath(args.input)
    barcode_design = os.path.normpath(args.barcodes)
    # output
    output_dir = os.path.normpath(args.output)
    exp_name = os.path.basename(output_dir)

    # validate user options 
    parser.validateArgs(args)

    if args.plasmid == 'a31': # (default option)
        # set plasmid path to default a.31 plasmid
        args.plasmid = str(files('barcodezy.plasmids').joinpath('a.31.fa'))
        # set flanks path to default a.31 flanks
        args.flanks = str(files('barcodezy.plasmids').joinpath('a.31_flanks.fa'))

    # remove any adaptor sequences (requires porechop), condense reads into a single file
    trimmed_read_file = f'{output_dir}/{exp_name}_trimmed.fastq'
    if not os.path.exists(trimmed_read_file):
        cmd_str = (f'porechop -i {read_dir} -o {trimmed_read_file} ' 
                   '--no_split --min_trim_size 25 '
                   '--adapter_threshold 97 --check_reads 2000')
        os.system(cmd_str)
        preprocess.rename_reads(trimmed_read_file, exp_name)

    # generate ref if SBARRO (also save ref_len for html plotting)
    if args.SBARRO:
        ref_len = preprocess.generate_ref(output_dir, args.insert_length)
    # else open ref fasta and concatenate to mimic circularity
    else:
        ref_len = preprocess.process_ref(output_dir, args.plasmid)

    # map trimmed reads to plasmid with minimap2
    summary_align_counts = preprocess.mm2_align(output_dir, trimmed_read_file, args.a31)

    # generate barcode alignment & read stats written to csv output
    align = barcode_aligner.barcode_scores(output_dir, barcode_design, args.flanks,
                                           args.insert_length, args.enzymes, 
                                           args.a31, args.SBARRO)
    # output compressed alignment csv with full set of BC scores and other metrics 
    align.to_csv(f'{output_dir}/{exp_name}.csv.gz', index=False)

    # analysis html and processing of alignment table (z-scores, top BC call, etc) 
    report = analysis.report_gen(output_dir, align, summary_align_counts, ref_len, args.zscore)
    # output compressed alignment csv with full set of BC scores and other metrics 
    report.to_csv(f'{output_dir}/{exp_name}_summary.csv.gz')

if __name__ == "__main__":
    main()