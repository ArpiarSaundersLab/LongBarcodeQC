import argparse
import sys
import os

DEFAULT_PLASMID = 'a31'

def getArgs() -> tuple[argparse.Namespace, argparse.ArgumentParser]:
    """Parse and return command-line arguments for barcodezy.

    Returns a tuple of (argparse.Namespace, ArgumentParser) for use with validateArgs.
    """
    arg_parser = argparse.ArgumentParser(
        description="Analyze long-read barcodes and generate a summary report"
    )
    arg_parser.add_argument('-i', '--input',
                        help='Path to input fastq files '
                        '(nanopore fastq_pass folder). ',
                        type=str,
                        required=True)
    arg_parser.add_argument('-o', '--output',
                        help='Path to output directory. '
                        'An output folder will be created if it does not '
                        'exist. Program will exit if files already exist here.',
                        type=str,
                        required=True)
    arg_parser.add_argument('-p', '--plasmid',
                        help='Path to plasmid fasta file. '
                        'The default plasmid is a.31. If using another '
                        'plasmid, use this option with the path to a fasta '
                        'file containing your plasmid. For best results, '
                        'concatenate the circular plasmid sequence to itself '
                        'to optimize alignment.',
                        type=str,
                        default=DEFAULT_PLASMID)
    arg_parser.add_argument('-b', '--barcodes',
                        help='Barcode FASTA file. Use a preset keyword or provide a path to a '
                        'custom FASTA file. Preset options: EV (Expression Vector), '
                        'AP (Assembly Plasmid), TS (TritSeq). '
                        'Custom file: path to a FASTA file with one barcode per entry; '
                        'headers must be unique.',
                        type=str,
                        required=True)
    arg_parser.add_argument('-l', '--insert_length',
                        help='Integer length of expected insert size cloned into MCS. '
                        'A rough estimate is sufficient. Default: 300 bp',
                        type=int,
                        default=300)
    arg_parser.add_argument('-f', '--flanks',
                        help='Path to fasta file containing MCS flanking regions. '
                        'The first sequence must be the upstream flank, and the second '
                        'sequence must be the downstream flank. '
                        'Must be provided in the same file.',
                        type=str)
    arg_parser.add_argument('-r', '--enzymes',
                        help='Path to text file with restriction site names and their '
                        'sequences. Define one site per line, with the name and sequence '
                        'separated by comma.',
                        type=str)
    arg_parser.add_argument('-a', '--a31',
                        help='Option to add a31 alignment. Mainly used '
                        'to detect contamination from a31 plasmid. This will flag a31 reads '
                        'in the final output and score barcodes against them.',
                        action='store_true')
    arg_parser.add_argument('-S', '--SBARRO',
                        help='Use this option if the plasmid is part of the SBARRO system '
                        '(derived from c.18). Reference will be generated with '
                        'NNN sequence inserted into the MCS (length of NNNs equal to '
                        'provided insert length).',
                        action='store_true')
    arg_parser.add_argument('-z', '--zscore',
                        help='Optionally set z-score threshold for barcode calling. The '
                        'default threshold is set by using the z-score for the worst '
                        'performing barcode alignment in the library (absolute value). '
                        'Sometimes this value is either too conservative or too permissive. '
                        'See the z-score distribution in the html output summary to evaluate this. '
                        'If desired - use this option to manually set the z-score threshold to be a '
                        'fixed value (e.g. 2.5).',
                        type=float)
    arg_parser.add_argument('-N', '--expected_insertions',
                        help='Optionally set the expected number of insertions in a library. '
                        'This is used in the html report to categorize reads in the read length histogram. '
                        'By default, this is set to the maximum number of sites/positions found '
                        'in the dataset.',
                        type=int)
    arg_parser.add_argument('--full-output',
                        help='Write the full barcode alignment scores as a parquet file '
                        '(includes MCS sequences and all per-barcode scores). '
                        'This file can be large.',
                        action='store_true')
    return arg_parser.parse_args(), arg_parser

def validateArgs(args: argparse.Namespace, arg_parser: argparse.ArgumentParser) -> None:
    """Validate arguments and filesystem state.

    - Ensures input/output paths exist in expected forms
    - Validates combinations of `--plasmid`, `--a31`, and `--SBARRO`
    - Checks barcode and restriction enzyme files if provided
    """
    input_path = os.path.normpath(args.input)
    output_path = os.path.normpath(args.output)
    barcode_path = os.path.normpath(args.barcodes)

    # exit if input path does not exist
    if not os.path.isdir(input_path):
        arg_parser.error(f'Input path does not exist or is not a directory: {input_path}')

    # create output directory if it does not exist
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # exit if output directory already contains files
    if len(os.listdir(output_path)) > 0:
        arg_parser.error('Output path is not empty. Remove contents or choose a different output name.')

    # verify reference options
    if args.plasmid != DEFAULT_PLASMID:  # if not a31, user provided a plasmid path
        if not os.path.exists(os.path.normpath(args.plasmid)):
            arg_parser.error(f'Plasmid fasta file does not exist: {args.plasmid}')
        if args.SBARRO:  # don't combine SBARRO and custom plasmid options
            arg_parser.error('SBARRO option cannot be used with a custom plasmid.')
    # cannot combine default a.31 plasmid and a.31 option (unless using SBARRO option)
    elif args.a31 and not args.SBARRO:
        arg_parser.error('No plasmid was provided. The default a.31 plasmid cannot be combined with '
                         'the a.31 option. Only use the a.31 option to annotate it as contamination.')

    BARCODE_PRESETS = {'EV', 'AP', 'TS'}
    if args.barcodes not in BARCODE_PRESETS:
        if not os.path.exists(barcode_path):
            arg_parser.error(
                f'--barcodes must be a preset (EV, AP, TS) or a valid file path. '
                f'File not found: {barcode_path}'
            )

    # verify restriction enzyme path and format
    if args.enzymes:
        rs_path = os.path.normpath(args.enzymes)
        if not os.path.exists(rs_path):
            arg_parser.error(f'Restriction enzyme file does not exist: {rs_path}')
        with open(rs_path, 'r') as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) != 2:  # check if each line has 2 elements (name, seq)
                    arg_parser.error('Restriction enzyme file must contain one name and sequence '
                                     'per line, separated by a comma.')
                seq = parts[1]
                if not all(c in {'A', 'C', 'T', 'G'} for c in seq.upper()):  # verify seq is only ACTG
                    arg_parser.error('Provided restriction enzyme sequences must only contain [ACTG]')

    # flanks required when using a custom plasmid without SBARRO
    if not args.SBARRO and args.plasmid != DEFAULT_PLASMID:
        if not args.flanks:
            arg_parser.error('Upstream and downstream MCS flanking regions must be provided (fasta), '
                             'unless using default a.31 plasmid or SBARRO option.')
        if not os.path.exists(os.path.normpath(args.flanks)):
            arg_parser.error(f'Flanks fasta file does not exist: {args.flanks}')

    # print user options
    print('User command:\n', ' '.join(sys.argv), '\n', sep='')
