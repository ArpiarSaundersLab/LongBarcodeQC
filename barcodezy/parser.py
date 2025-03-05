import argparse
import sys, os

def getArgs():
    parser = argparse.ArgumentParser()
    # arg options
    parser.add_argument('-i', '--input',
                        help='Path to input fastq files '
                        '(nanopore fastq_pass folder). ',
                        type=str,
                        required=True)
    parser.add_argument('-o', '--output',
                        help='Path to output directory. '
                        'An output folder will be created if it does not '
                        'exist. Program will exit if files already exist here.',
                        type=str,
                        required=True)
    parser.add_argument('-p', '--plasmid',
                        help='Path to plasmid fasta file. '
                        'The default plasmid is a.31. If using another '
                        'plasmid, use this option with the path to a fasta '
                        'file containing your plasmid. For best results, '
                        'concatenate the circular plasmid sequence to itself '
                        'to optimize alignment.',
                        type=str,
                        default='a31')
    parser.add_argument('-b', '--barcodes',
                        help='Path to barcode fasta file. '
                        'One barcode per line. Fasta headers must be '
                        'unique names.',
                        type=str,
                        required=True)
    parser.add_argument('-l', '--insert_length',
                        help='Integer length of expected insert size cloned into MCS. '
                        'A rough estimate is sufficient. Default: 300 bp',
                        type=int,
                        default=300)
    parser.add_argument('-f', '--flanks',
                        help='Path to fasta file containing MCS flanking regions. '
                        'The first sequence must be the upstream flank, and the second '
                        'sequence must be the downstream flank. '
                        'Must be provided in the same file.',
                        type=str)
    parser.add_argument('-r', '--enzymes',
                        help='Path to text file with restriction site names and their '
                        'sequences. Define one site per line, with the name and sequence '
                        'separated by comma.',
                        type=str)
    parser.add_argument('-a', '--a31',
                        help='Option to add a31 alignment. Mainly used '
                        'to detect contamination from a31 plasmid. This will flag a31 reads '
                        'in the final output and score barcodes against them.',
                        action='store_true')
    parser.add_argument('-S', '--SBARRO',
                        help='Use this option if the plasmid is part of the SBARRO system '
                        '(derived from c.18). Reference will be generated with '
                        'NNN sequence inserted into the MCS (length of NNNs equal to '
                        'provided insert length).',
                        action='store_true')
    parser.add_argument('-z', '--zscore',
                        help='Optionally set z-score threshold for barcode calling. The '
                        'default threshold is set by using the z-score for the worst '
                        'performing barcode alignment in the library (absolute value). '
                        'Sometimes this value is either too conservative or too permissive. '
                        'See the z-score distribution in the html output summary to evaluate this. '
                        'If desired - use this option to manually set the z-score threshold to be a '
                        'fixed value (e.g. 2.5).',
                        type=float)
    return parser.parse_args()

def validateArgs(args) -> None:
    input_path = os.path.normpath(args.input)
    output_path = os.path.normpath(args.output)
    plasmid_path = args.plasmid
    barcode_path = os.path.normpath(args.barcodes)
    restriction_path = args.enzymes

    # exit if input path does not exist
    if not os.path.isdir(input_path):
        print(f'Error: Input path does not exist or is not a directory: {input_path}')
        sys.exit(1)

    # create output directory if it does not exist
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # exit if output directory already contains files
    if len(os.listdir(output_path)) > 0:
        # allow <= 1 file to exist to avoid regenerating the trimmed fastq output
        print('Error: Output path is not empty. Remove contents or choose a different output name.')
        sys.exit(1)
    
    # verify reference options
    if args.plasmid != 'a31': # if not a31, user input plasmid path 
        if not os.path.exists(os.path.normpath(plasmid_path)):
            print(f'Error: Plasmid fasta file does not exist: {plasmid_path}')
            sys.exit(1)
        if args.SBARRO: # don't combine SBARRO and plasmid reference options
            print('Error: SBARRO option cannot be used with a custom plasmid.')
            sys.exit(1)
    # cannot combine default a.31 plasmid and a.31 option (unless using SBARRO option)
    elif args.a31 and not args.SBARRO: 
        print('Error: No plasmid was provided. The default a.31 plasmid cannot be combined with '
              'the a.31 option. Only use the a.31 option to annotate it as contamination.')
        sys.exit(1)

    if not os.path.exists(barcode_path): # check if barcode path exists
        print(f'Error: Barcode fasta file does not exist: {barcode_path}')
        sys.exit(1)
    
    # verify restriction enzyme path exists
    if args.enzymes:
        rs_path = os.path.normpath(restriction_path)
        if not os.path.exists(rs_path):
            print(f'Error: Restriction enzyme file does not exist: {rs_path}')
            sys.exit(1)
        with open(rs_path, 'r') as fh:
            for rs in fh:
                rs = rs.strip().split(',') 
                if len(rs) != 2: # check if each line has 2 elements (name, seq)
                    print('Error: Restriction enzyme file must contain one name and sequence '
                          'per line, separated by a comma.')
                    sys.exit(1)

                seq = rs.strip().split(',')[1]
                # very seq is only ACTG
                if not all(c in {'A', 'C', 'T', 'G'} for c in seq):
                    print('Error: Provided restriction enzyme sequences must only contain [ACTG]')
                    sys.exit(1)
    
    # plasmid path and flanks must be provided without SBARRO option
    if not args.SBARRO and args.plasmid != 'a31':
        # require flanking sequence around MCS if not default a.31 or SBARRO
        if not args.flanks:
            print('Error: Upstream and downstream MCS flanking regions must be provided (fasta), '
                  'unless using default a.31 plasmid or SBARRO option.')
            sys.exit(1)