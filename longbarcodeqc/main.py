#!/usr/bin/env python
import glob
import os
import sys
from datetime import datetime
from importlib.resources import files
from longbarcodeqc import analysis
from longbarcodeqc import barcode_aligner
from longbarcodeqc import consensus
from longbarcodeqc import parser
from longbarcodeqc import preprocess
from longbarcodeqc import trim


def main(args=None):
    """Entry point for the LongBarcodeQC CLI workflow."""
    start = datetime.now()

    args, arg_parser = parser.getArgs()
    print(f'\nStarting run...\n{start.strftime("%Y-%m-%d %H:%M:%S")}\n')

    read_dir = os.path.normpath(args.input)
    barcode_design = os.path.normpath(args.barcodes)
    output_dir = os.path.normpath(args.output)
    exp_name = os.path.basename(output_dir)

    # validate user options (this creates the output directory)
    parser.validateArgs(args, arg_parser)

    # record the run details at the top of the log
    user_command = ' '.join(sys.argv)
    with open(f'{output_dir}/Log.txt', 'a') as log:
        log.write(f'LongBarcodeQC version: {parser.getVersion()}\n')
        log.write(f'Run started: {start.strftime("%Y-%m-%d %H:%M:%S")}\n')
        log.write(f'User command: {user_command}\n\n')

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
    if args.SBARRO:
        args.flanks = str(files('longbarcodeqc.plasmids').joinpath('SBARRO_flanks.fa'))

    read_file = f'{output_dir}/{exp_name}.fastq.gz'
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
    # keep the target plasmid sequence for a -p consensus reference
    # (the concatenated ref file is removed after alignment)
    ref_seq = preprocess.read_ref(output_dir) if args.consensus else None

    # map reads to plasmid with minimap2
    summary_align_counts = preprocess.mm2_align(output_dir, read_file, args.AP, is_default_plasmid,
                                                reads_removed_in_trimming)

    # generate barcode alignment & read stats written to csv output
    backbone_prefix = f'{output_dir}/.{exp_name}.backbone_reads' if args.consensus else None
    align = barcode_aligner.barcode_scores(output_dir, barcode_design, args.flanks,
                                           args.insert_length, args.enzymes,
                                           args.AP, is_default_plasmid, args.SBARRO,
                                           backbone_prefix)

    # optional consensus of the plasmid backbone (everything outside the MCS, flanks included),
    # compared to the matching reference backbone:
    #   -S: the expression vector (EV) map, whatever the barcode set
    #   default AP plasmids: AP-Amp or AP-Kan, whichever has more reads (built from those reads)
    #   -p: the provided plasmid, if it contains both flanks
    backbone = None
    if args.consensus:
        if args.SBARRO:
            backbone_type, reference = 'plasmid', consensus.load_reference('EV')
        elif is_default_plasmid:
            backbone_type = ('AP-Kan' if summary_align_counts.get('AP-Kan', 0) >
                             summary_align_counts.get('AP-Amp', 0) else 'AP-Amp')
            reference = consensus.load_reference(backbone_type)
        else:
            backbone_type = 'plasmid'
            reference = consensus.reference_from_plasmid(ref_seq, args.flanks, args.insert_length,
                                                         os.path.basename(args.plasmid))
        print(f'Building backbone consensus ({reference["name"] if reference else backbone_type})...')
        backbone = consensus.backbone_consensus(output_dir, f'{backbone_prefix}.{backbone_type}.fa',
                                                args.flanks, args.insert_length, reference)
        if backbone is not None:
            backbone['read_type'] = backbone_type
        for leftover in glob.glob(f'{backbone_prefix}.*.fa'):
            os.remove(leftover)

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
        user_command,
        backbone,
        parser.getVersion(),
    )
    print('Writing summary csv...')
    report.drop(columns=['MCS_seq']).to_csv(f'{output_dir}/{exp_name}_summary.csv.gz')

    print(f'\nRun complete!')
    end = datetime.now()
    print(end.strftime("%Y-%m-%d %H:%M:%S"))
    print(f'Time elapsed: {end - start}')

    with open(f'{output_dir}/Log.txt', 'a') as log:
        log.write(f'\nRun finished: {end.strftime("%Y-%m-%d %H:%M:%S")}\n')
        log.write(f'Time elapsed: {end - start}\n')


if __name__ == "__main__":
    main()
