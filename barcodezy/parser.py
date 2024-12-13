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
                        required=False)
    parser.add_argument('-b', '--barcodes',
                        help='Path to barcode fasta file. '
                        'One barcode per line. Fasta headers must be '
                        'unique names.',
                        type=str,
                        required=True)
    parser.add_argument('-l', '--insert_length',
                        help='Integer length of expected insert size cloned into MCS.',
                        type=int,
                        required=True)
    parser.add_argument('-f', '--flanks',
                        help='Path to fasta file containing MCS flanking regions. '
                        'The first sequence must be the upstream flank, and the second '
                        'sequence must be the downstream flank. '
                        'Must be provided in the same file.',
                        type=str,
                        required=False)
    parser.add_argument('-r', '--enzymes',
                        help='Path to text file with restriction site names and their '
                        'sequences. Define one site per line, with the name and sequence '
                        'separated by comma.',
                        type=str,
                        required=False)
    parser.add_argument('-a', '--a31',
                        help='Option to add a31 alignment. Mainly used '
                        'to detect contamination from a31 plasmid. This will flag a31 reads '
                        'in the final output and score barcodes against them.',
                        action='store_true',
                        required=False)
    parser.add_argument('-S', '--SBARRO',
                        help='Use this option if the plasmid is part of the SBARRO system '
                        '(derived from c.18). Reference will be generated with '
                        'NNN sequence inserted into the MCS (length of NNNs equal to '
                        'provided insert length).',
                        action='store_true',
                        required=False)
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
    if len(os.listdir(output_path)) > 1:
        # allow <= 1 file to exist to avoid regenerating the trimmed fastq output
        print('Error: Output path is not empty. Remove contents or choose a different output name.')
        sys.exit(1)
    
    # verify plasmid and barcode fasta files exist
    if args.plasmid:
        if not os.path.exists(os.path.normpath(plasmid_path)):
            print(f'Error: Plasmid fasta file does not exist: {plasmid_path}')
            sys.exit(1)
    if not os.path.exists(barcode_path):
        print(f'Error: Barcode fasta file does not exist: {barcode_path}')
        sys.exit(1)
    
    # verify restriction enzyme path exists
    if args.enzymes:
        if not os.path.exists(os.path.normpath(restriction_path)):
            print(f'Error: Restriction enzyme file does not exist: {restriction_path}')
            sys.exit(1)
    
    # plasmid path and flanks must be provided without SBARRO option
    if not args.SBARRO:
        # require plasmid reference fasta file
        if not args.plasmid:
            print('Error: Plasmid reference fasta file must be provided. Otherwise, '
                  'use the SBARRO option to generate the default reference.')
            sys.exit(1)
        # require flanking sequence around MCS
        if not args.flanks:
            print('Error: Upstream and downstream MCS flanking regions must be provided (fasta).')
            sys.exit(1)