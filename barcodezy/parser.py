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
                        required=False,
                        default='a.31')
    parser.add_argument('-b', '--barcodes',
                        help='Path to barcode fasta file. '
                        'One barcode per line. Fasta headers must be '
                        'unique names.',
                        type=str,
                        required=True)
    parser.add_argument('-l', '--insert_length',
                        help='Integer length of expected insert size cloned into MCS ',
                        type=int,
                        required=True)
    return parser.parse_args()

def validateArgs(input_path, output_path, plasmid_path, barcode_path) -> None:
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
    if not os.path.exists(plasmid_path):
        print(f'Error:Plasmid fasta file does not exist: {plasmid_path}')
        sys.exit(1)
    if not os.path.exists(barcode_path):
        print(f'Error: Barcode fasta file does not exist: {barcode_path}')
        sys.exit(1)