#!/usr/bin/env python
import os
from datetime import datetime
from importlib.resources import files
from longbarcodeqc import analysis
from longbarcodeqc import barcode_aligner
from longbarcodeqc import parser
from longbarcodeqc import preprocess
from longbarcodeqc import trim


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

    BARCODE_PRESETS = {
        'EV': 'barcodes/3_1_256_rev/barcodes_tail_only.fa',
        'AP': 'barcodes/3_1_256/barcodes_tail_only.fa',
        'TS': 'barcodes/4_4_3/barcodes_probe_only.fa',
    }
    if args.barcodes in BARCODE_PRESETS:
        barcode_design = str(files('longbarcodeqc').joinpath(BARCODE_PRESETS[args.barcodes]))

    is_default_plasmid = (args.plasmid == parser.DEFAULT_PLASMID) and not args.SBARRO
    if args.plasmid == parser.DEFAULT_PLASMID: # (default option)
        # set plasmid path to default AP-Amp plasmid
        args.plasmid = str(files('longbarcodeqc.plasmids').joinpath('AP-Amp.fa'))
        # set flanks path to default AP flanks
        args.flanks = str(files('longbarcodeqc.plasmids').joinpath('AP_flanks.fa'))

    read_file = f'{output_dir}/{exp_name}.fastq'
    preprocess.rename_reads(read_dir, read_file, exp_name)

    # trim ONT Rapid adapter before alignment
    reads_removed_in_trimming = None
    if args.trim:
        print('Trimming adapters...')
        reads_removed_in_trimming = trim.trim_fastq_in_place(
            read_file,
            log_file=f'{output_dir}/{exp_name}.cutadapt.txt',
            threads=3,
        )

    if args.SBARRO:
        preprocess.generate_ref(output_dir, args.insert_length)
    else:
        preprocess.process_ref(output_dir, args.plasmid)
    print('Preparing reference plasmid...')

    # map reads to plasmid with minimap2
    summary_align_counts = preprocess.mm2_align(output_dir, read_file, args.AP, is_default_plasmid,
                                                reads_removed_in_trimming)

    # generate barcode alignment & read stats written to csv output
    align = barcode_aligner.barcode_scores(output_dir, barcode_design, args.flanks,
                                           args.insert_length, args.enzymes,
                                           args.AP, is_default_plasmid, args.SBARRO)
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
