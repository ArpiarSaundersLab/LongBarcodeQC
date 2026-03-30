#!/usr/bin/env python
import os
from datetime import datetime
from importlib.resources import files
from longbarcodeqc import analysis
from longbarcodeqc import barcode_aligner
from longbarcodeqc import parser
from longbarcodeqc import preprocess


def main(args=None):
    """Entry point for the LongBarcodeQC CLI workflow."""
    start = datetime.now()
    print(f'\nStarting run...\n{start.strftime("%Y-%m-%d %H:%M:%S")}\n')

    args, arg_parser = parser.getArgs()

    read_dir = os.path.normpath(args.input)
    barcode_design = os.path.normpath(args.barcodes)
    output_dir = os.path.normpath(args.output)
    exp_name = os.path.basename(output_dir)

    # validate user options
    parser.validateArgs(args, arg_parser)

    if args.plasmid == parser.DEFAULT_PLASMID: # (default option)
        # set plasmid path to default a.31 plasmid
        args.plasmid = str(files('longbarcodeqc.plasmids').joinpath('a.31.fa'))
        # set flanks path to default a.31 flanks
        args.flanks = str(files('longbarcodeqc.plasmids').joinpath('a.31_flanks.fa'))

    read_file = f'{output_dir}/{exp_name}.fastq'
    preprocess.rename_reads(read_dir, read_file, exp_name)

    # generate ref if SBARRO (also save ref_len for html plotting)
    if args.SBARRO:
        ref_len = preprocess.generate_ref(output_dir, args.insert_length)
    # else open ref fasta and concatenate to mimic circularity
    else:
        ref_len = preprocess.process_ref(output_dir, args.plasmid, args.insert_length)
    print(f'Preparing reference plasmid ({ref_len} bp)...')

    # map reads to plasmid with minimap2
    summary_align_counts = preprocess.mm2_align(output_dir, read_file, args.a31)

    # generate barcode alignment & read stats written to csv output
    align = barcode_aligner.barcode_scores(output_dir, barcode_design, args.flanks,
                                           args.insert_length, args.enzymes, 
                                           args.a31, args.SBARRO)
    # output verbose parquet with every barcode alignment score per read (can be large)
    if args.full_output:
        print('\nWriting full barcode alignment file...')
        align.to_parquet(f'{output_dir}/{exp_name}.parquet', index=False)

    # analysis html and processing of alignment table (z-scores, top BC call, etc)
    print('Generating html report...')
    report = analysis.report_gen(
        output_dir,
        align,
        summary_align_counts,
        ref_len,
        args.zscore,
        args.expected_insertions,
    )
    print('Writing summary csv...')
    report.drop(columns=['MCS_seq']).to_csv(f'{output_dir}/{exp_name}_summary.csv.gz')

    print(f'\nRun complete!')
    end = datetime.now()
    print(end.strftime("%Y-%m-%d %H:%M:%S"))
    print(f'Time elapsed: {end - start}')


if __name__ == "__main__":
    main()
