import glob
import gzip
import os
from importlib.resources import files
from typing import Dict


def rename_reads(input_path: str, output_file: str, exp_name: str) -> None:
    """Merge FASTQ files from a directory (or single file) and rewrite headers
    with experiment name and sequential IDs."""
    if os.path.isfile(input_path):
        fastq_files = [input_path]
    else:
        fastq_files = sorted(
            glob.glob(os.path.join(input_path, '**', '*.fastq'), recursive=True) +
            glob.glob(os.path.join(input_path, '**', '*.fastq.gz'), recursive=True) +
            glob.glob(os.path.join(input_path, '**', '*.fq'), recursive=True) +
            glob.glob(os.path.join(input_path, '**', '*.fq.gz'), recursive=True)
        )
        if not fastq_files:
            raise FileNotFoundError(f'No FASTQ files found in: {input_path}')

    print(f'Merging and renaming reads from {len(fastq_files)} file(s)...\n')
    read_count = 0
    with open(output_file, 'w') as out:
        for fq_path in fastq_files:
            opener = gzip.open if fq_path.endswith('.gz') else open
            with opener(fq_path, 'rt') as reads:
                for i, line in enumerate(reads):
                    if i % 4 == 0:
                        read_count += 1
                        parts = line.strip().split(' ')
                        rest = ' '.join(parts[1:])
                        out.write(f'@{exp_name}_{read_count} {rest}\n')
                    else:
                        out.write(f'{line.strip()}\n')

def generate_ref(outpath: str, insert_len: int) -> int:
    exp_name = os.path.basename(outpath)
    upstream_MCS = 'acgcttaacaacaaaaccagagaagaaaaagacagcgtcaattgcaaagcaaaaatgtaacacccctacaatggatgccgacaagattgtgttcaaagtcaataatcaggtggtctctttgaagcctgagattatcgtggatcaatatgagtacaagtaccctgccatcaaggatttgaaaaagccttgtatcaccctagggaaagcccccgacttgaacaaagcatacaaatcagttttatcaggcatgaatgccgccaaacttgatccggatgatgtatgctcctacttggcagcagcaatgcagttctttgaggggacatgtccggaagactggaccagctatggaatcctgattgcacgaaaaggagataggatcaccccaaactctctagtggagataaagcgtactgatgtagaagggaattgggctctgacaggaggcatggaattgacaagggaccccactgtctctgaacatgcatctttagtcggtcttctcctgagtctgtacaggttgagcaaaatatcaggacagaacactggtaactataagacaaacattgcagataggatagagcagattttcgagacagcaccttttgttaagatcgtggaacaccataccctaatgacaactcacaagatgtgtgctaattggagtactataccgaacttcagatttttggccggaacctacgacatgtttttctcacggattgagcatctatattcggcaatcagagtgggcacagtcgtcaccgcttatgaagactgctcaggactggtatcgtttacagggttcataaagcagatcaatctcaccgcaagggaagcaatactatatttcttccacaagaactttgaggaagagataagaagaatgttcgagccagggcaagagacagctgttcctcactcttatttcatccacttccgttcactaggcttgagtgggaagtctccttattcatcgaatgctgtcggtcatgtgttcaatctcattcactttgttggatgctacatgggtcaagtcagatctctaaatgcgacggttattgctgcatgtgcccctcatgagatgtctgttctagggggctatttgggagaggaattcttcggaaaagggacatttgaaagaaggttcttcagagacgagaaagaacttcaagaatatgaggcggctgaactaacaaagtccgacgtggcactggcagatgacggaaccgtcaactctgatgacgaggactatttctctggtgaaaccagaagtccagaagctgtctatactcgaatcatgatgaatggaggtcgactgaagagatctcatatacggagatatgtctcagtcagttccaatcatcaagcccgtccaaactcattcgccgaatttttaaacaagacgtattcgaatgactcataaggagttgattgacagggtgccagaaatctatagattgtatatatccatcatgaaaaaaactaacactcctcctttcaaaccatcccaaatatgagcaagatctttgttaatccgagtgcaatcagagccggtctggccgatcttgagatggccgaagagactgttgatctgatcaacagaaacatagaagacaatcaggctcatctccagggagaacccatagaagtggacaacctccctgaggacatgaagcgacttcacctggacgatgaaaaatcgtccaaccttggtgagatggttagggtgggagaaggcaagtatcgagaggactttcagatggatgagggagaggaccccaacctcctgttccaatcgtacctggataatgttggagtccaaatagtcagacaaatgaggtcaggagagagattcctcaagatatggtcacagaccgtagaggaaattgtatcctatgtcacggtcaactttcctaaccctccaagaaggtcttcggaggataaatcaacccagactactggcagagagctcaagaaggagacaacgtctgctttctctcagagagaaagccaaccttcgaaagctaggatggtggctcaagttgcccctggtcctccagcccttgaatggtcagccaccaatgaagaagatgatctatcagtagaggctgagatcgctcatcagattgctgaaagcttttccaagaagtacaagtttccctcccgatcttcaggaatattcttgtataattttgagcaactgaagatgaaccttgatgacatagttaaagaggcaaaaaatgtaccgggcgtgacccgtctggcccatgatggatccaaaatccccctgagatgcgtactgggatgggtcgctttggctaattccaaaaaattccaattactagtcgaggctgacaagctaagcaaaatcatgcaagatgatttgaatcgctatacatcctgctaaccgagttttcgaactcagtccctccagataatgaaaactgagatgttatggagtagacatgaaaaaaacaggcaacaccactgataaaatgaacgttctacgcaagatagtgaaaaaatgtagggatgaggacactcaaaagccctctcctgtgtcagcccctccgtatgacgatgacctgtggcttccacctcctgaatatgtcccgctgaaagaactcacaagcaagaagaacatgaggaacttttgtgtcaacggggaggttaaagcgtgtagcccaaatggttactcattcaggattttgcggcacattctgagatcattcaacgagatatactctgggaatcataggatgattgggttagtcaaaattgttgttggactagctttatcaggagctccagtacctgagggcatgaactgggtatacaaattgaggagaacccttatattccagtgggctgattccaggggccctcttgaaggggaggagttagattactctcaagagatcacttgggatgatgatactgaattcgtcggattgcaaataagagtgagcgcaagacaatgtcatattcaaggcaggatctggtgtatcaacacgaactcgagggcatgtcaactatggtctgacatgtctcttcagacacaaaggtctgaagaggacaaagactcttctctgcttctagaataatcagattatatcctgcaagtgtatcacttgtttacctctggaggagagagcatacaggcttgactccgatccttgggagcaatagaacaaaaaaacacacgttatggtgccgttaaatcgctgcattttatcaaagtcaagttgataacctttacattttgagcctcttggatgtgaaaaaaactattaacatccctcaaaagaccccggtaaagccaccatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaaggctcgagctgatccaaaaaagaagagaaaggtagatccaaaaaagaagagaaaggtagatccaaaaaagaagagaaaggtaggatccagataactgatcataatcagccataccacatttggaattcgtcgagggggc'.upper()
    downstream_MCS = 'gctagacatgaaaaaaactaacactcctccggtaccgccaccatggacgtggtgaatcagctggtggctgggggtcagttccgggtggtcaaggagccccttggcttcgtgaaggtgctgcagtgggtctttgccatcttcgcctttgctacgtgtggcagctacaccggggagcttcggctgagcgtggagtgtgccaacaagacggagagtgccctcaacatcgaagttgaattcgagtaccccttcaggctgcaccaagtgtactttgatgcaccctcctgcgtcaaagggggcactaccaagatcttcctggttggggactactcctcgtcggctgaattctttgtcaccgtggctgtgtttgccttcctctactccatgggggccctggccacctacatcttcctgcagaacaagtaccgagagaacaacaaagggcctatgatggactttctggctacagccgtgttcgctttcatgtggctagttagttcatcagcctgggccaaaggcctgtccgatgtgaagatggccacggacccagagaacattatcaaggagatgcccatgtgccgccagacagggaacacatgcaaggaactgagggaccctgtgacttcaggactcaacacctcagtggtgtttggcttcctgaacctggtgctctgggttggcaacttatggttcgtgttcaaggagacaggctgggcagccccattcatgcgcgcacctccaggcgccccggaaaagcaaccagcacctggcgatgcctacggcgatgcgggctacgggcagggccccggaggctatgggccccaggactcctacgggcctcagggtggttatcaacccgattacgggcagccagccagcggtggcggtggctacgggcctcagggcgactatgggcagcaaggctatggccaacagggtgcgcccacctccttctccaatcagatgggatccatcgccaccatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaagctagcggcgcctagccggtcatccttttgacgattccagtcccgaggataacctcctctcgggattggggggaatctttggatccagtagtcctccttgaactccgtccaacagggcagattcaagagtcataagactttcattaatcatttcagttgatcagacatggtcgtgtagattctcttaatacgggagatcttctagcagtttcagtgaccaacggtgctttcattctccaggaactgataccaaaggttgtggacaggccaaggggtacttcggatgactctgtgcttgggcacagaaagaggtcgtagtgtgccccccgatagcggactcgacatgaatcaactaagaaaggcaatctgcctcccatgatggacataagcaatagttcacaaccatcttgcatctcagtgaagtgtacataactattgagggctgggtcatctaagcatttcagtcgagaaaaaaactgtagaccaaaagaacaactagcaacacttctcatccagagacccatatcaagatgctagatccgggagaggtttatgatgaccctattgatccaattgagtcagaggctgaacccagaggaacccccactgtccccaacatcttgaggaactccgactacaatctcaattctcctttgatagaggatcctgccaaactaatgttagaatggttgaagacagggaacagaccttatcggatgactttgacagacaattgctccaggtcttacaaagttttgaaagattatttcaagaaagtagatttgggttctctcaaagtgggcggaactgctgcacagtcaatggtttctctctggttgtgtggtgcccactctgaatcaaacaggagccggagatgtataaccgacttggcccatttctattccaagtcatcccccatagagaagctattgaattgtacgctaggaaacagaggcctgagaatcccaccagagggggtgttaaattgcctcgagagggtcaattatgacaaggcatttgggaggtatctggccaacacgtattcctcttacttgtttttccatgtaatcaccttatacatgaatgccttagactgggaagaggaaaagaccatcctggcattatggaaagatataacctcagtggataccgagaaggacttggtcaaattcaaagatcaaatatggggactgttgattgtgacaaaggactttgtttactctcagagttcaaactgtctttttgacagaaactacacactgatgctaaaggatcttttcttgtctcgattcaactccttaatgattctgctttctccccctgagccccgatactcagatgacttaatatctcagctgtgccagctatacatcgctggggatcaagtcttgtccatgtgtgggaactccggctatgaagtcatcaaaatattggagccatatgtcgtgaacagtttggtccagagggcagagaagtttaggcctctcatccaccccttgggagactttcctatgttcataaaagacaaggtgaatcaacttgaagggacttttggtcccagtgcaaagaggttttttagggttttagatcaattcgacaacatacatgacttagtatttgtgtatggctgttacagacattgggggcacccctatatagattatcggaagggtctgtcgaaactatatgatcaagttcacattaagaaagtaatagataagtcctaccaggagtgtttagcaagtgacttggccagaaggatcctcagatggggatttgacaagtactccaaatggtatctagattcgagattccttgcccgagaccaccccttgactccttatgtcaagacccaaacatggccacccaaacatatagtagacttggtgggggacacatggcataagctcccgatcacgcagatctttgaaattcctgaatcaatggacccgtcagagatactggatgacaaatcacattctttcaccagaacaaggttagcttcttggctgtcagagaaccgaggggggcctgttcctagcgagaaggtcattatcacggccctgtctaagccacctgtcaatccccgagagtttttgaaatctatagacctcggaggattgccagatgaggatttgataattggcctcaaaccaaaggaacgggagttgaagatcgagggccgattctttgctctaatgtcatggaatctaagattatattttgtcatcaccgaaaagctcttggccaactacattttgccactttttgacgcactgactatgacagacaacctgaacaaggtgttcaaaaagttgatcgacagggtcaccgggcaagggcttttggactattctagggtcacatacgcatttcacctggactatgagaaatggaacaatcatcaaagattggagtcaacagaggatgtattctctgtcctagatcaggtgtttggattgaagagggtgttttctagaacacacgagttttttcagaagtcctggatctattattcagacagatcagacctcattgggttatgggaggaccaaatatattgcttggatatgtctaacggcccaacctgctggaatggccaagatggcgggctagagggcttgcggcagaagggctggagtctagtcagtttattaatgatagatagagaatctcaaaccaggaacacaagaaccaagatactagctcaaggagacaaccaggttctgtgtcctacatacatgttgtcaccgggattgtctcaagaggggcttctctatgagttagagagcatatcaaggaatgcactctcaatataccgagctatcgaggaaggagcatctaagctggggctgatcatcaagaaggaagagaccatgtgtagttatgactttctcatatatgggaagacccccttatttcgaggcaacatattggtgcctgaatccaaaagatgggcccgagtctcttgcatctctaacgaccaaatagtcaacctcgccaatataatgtcgacagtatccaccaatgcgctgacagtggcacaacactctcaatctctgatcaaacctatgagggattttctgctcatgtcagtacaggcagttttccactacctgttgtttagcccaatcctaaaaggcagagtttataagattctgagtgctgaaggggagagctttctcctagccatgtcgcggataatctacctagatccttctttgggaggggtgtctggaatgtctctcgggaggttccatatacgtcagttctcagaccctgtctctgaagggttgtcattctggagagagatctggttaggctctcatgagtcctggattcacgcgttgtgtcaggaggccgggaaccccgatcttggagagagaacactcgagagcttcactcgccttttagaagatcctactaccttaaatatcaaaggaggggccagccctaccattctactcaaggatgctatcagaaaggctctgtacgacgaggtggacaaggtggaaaattcagagtttcgagaggcaatcctgctgtccaagacccatagagataactttatactctttttaaaatctgttgagcctctgttccctcgatttctcagtgagctcttcagttcgtctttcttgggaataccggagtcaatcattggactgatacaaaactcccggacaataagaaggcagtttagaaagagtctctcaagaactttagaagagtccttctacaactcagagatccacgggattaatcggataacccagacacctcaaagggtcggaagggtgtggccttgctcttcagagagggcagatctacttagggagatctcttgggggaggaaagtggtaggcacgacagttcctcacccttccgagatgttggagttgtttcccaaatcctccatttcctgcacttgtggagcaacagggggaggcaatcctagagtctctgtatcagtactcccgtccttcgatcagtcatttttctcacggggccccctaaagggatacttgggctcgtccacctccatgtcaacccagctattccatgcatgggaaaaagtcactaatgttcatgtggtgaaaagggctctatcgttaaaagaatctataaactggttcatcaataggaattccaatttggctcaaactctaattggaaacatcatgtctctgacaggccctgatttccctctagaagaggcccctgttttcaaacggacagggtcagccttgcataggttcaagtctgccagatacagcgaaggagggtattcttctgtttgccctaaccttctctcccatatctctgttagtacagacactatgtctgatttgacccaaaacgggaagaactatgatttcatgtttcagccattgatgctttatgcgcaaacatggacatcggaactggtacagagggatacaagacttagagactccacgtttcactggcaccttcggtgcaacagatgtgtgaggcccattgatgatataacactggaaacttctcagatcttcgagttcccggatgtgtcaaaaaggatatccaggatggtttctggagctgtccctcaatttcagaagcttcctgatatccgtctaagaccaggtgattttgaatccctaagtggtagagaaaagtctcgccatatagggtcagctcaggggctcttatactcaatcttagtagcaattcacgactcaggatacaatgatgggaccatcttccctgtcaacatatacggcaaagtttcccccagagactatttgagagggcttgcaagagggatcttgatagggtcctcgatttgcttcttgacacgaatgacaaatattaacattaacagacctcttgaattgatctcaggggtaatttcctatattctcctgaggctggataatcatccctctctgtatataatgcttagagaaccgtctcttagaggagaaatattctctatccctcagaaaatccccgccgcttacccaaccactatgagagaaggcaacagatcgattttgtgttacctccaacacgtgctacgctatgagcgagaggcaatcacggcgtccccggagaatgactggctgtggatcttctcagacttcagaagtgtgaaaatgacgtacttgaccctcattacctaccagtctcacctcctactccagagggttgagagaaacttgtctaagagtatgagagctactctgcgacaaatgagttccttaatgaggcaagtgctgggtgggcacggagaagataccttggagtcagacgatgacattcaacgattactaaaagactctttgcgaaagacaaggtgggtggatcaagaggtgcgccatgcagctagaaccatgaatggagattacagccccgacaagaaagtatcccacaaggcaggatgttcagaatgggtctgctctgctcaacagattgccgtctccacctcagccaacccggcccctgtctcagagcttgacattagggccctctctaagaggtttcaaaaccccttgatctcgggcctgagagtggttcagtgggcaaccggtgcccattataagcttaagcctattctagatgatctaaatgttttcccatctctctgtcttgtaattggagacgggtcagggggaatatcaagggcagttctcaacatgtttccagattctaagcttgtgttcaacagcctattggaggtgaatgatctgatggcttccggaacacatccactgcctccttcagcaatcatgagtggaggagatgacatcatctccagagtgatagactttgactcaatatgggagaaaccgtccgacctgaggaactcggccacctggagatacttccagtcggttcaaaaacaggtcaacatgtcgtatgacctcattatttgtgatgcagaagttactgatattgcatctatcaaccggataactctgttgatgtctgatttcgcattgtctatagatggaccactttatctggtcttcaaaacttacgggactatgctagtaaacccggactataaagctattcaacatctgtcaagagcgttcccttcggtcacagggtttgtaacccaagtaacttcatccttttcttctgagctatacctccggttctctaaacgaggaaagtttttcagggacgccgagtacttgacctcttccacccttcgagagatgagccttgtgttgttcaattgcagcagccccaaaagtgagatgcagagagctcgttccttaaactatcaagacctggtaaggggatttcctgaagagatcatatcaaatccttacaacgaaatgatcataactctgattgacaatgatgtagagtccttcctagtccacaagatggtggatgatcttgagctacagaggggaactctgtctaaagtggctatcattatatccatcatgatcgttttttccaatagagtcttcaacatttccaaacctttgactgaccccttgttctaccccccatctgatcctaaaatcctgaggcacttcaacatatgttgcagtactatgatgtatctatctaccgctttaggcgacgtccctaacttcgcaagacttcatgacctgtataacagacctataacttgttacttcagaaaacaagttattcgagggaatatttatctatcttggagttggtccgatgataccccagtgttcaagagagtagcctgtaattctagcttgagtctgtcatctcactggatcaggttgatctacaagatagtgaagactaccagactcattggcagcataaaagacctatcaggagaggtagaacgacatcttcatgggtataacagatggatcaccctcgaggatatccgatctagatcatccctactagattacagttgtttgtaagccggatattaccgaaagcctgtgcatgctaaaattcttgtatgatgcatcttgaaaaaaacaagatcttgaatccggacctctggttgtttgattgttttttccatctttattgtttttttgttaagcgtgggtcggcatggcatctccacctcctcgcggtccgacctgggcatccgaaggaggacgcacgtccactcggatggctaagggagggcggggatccggctgctaacaaagcccgaaaggaagctgagttggctgctgccaccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgctgaaagtcgcgcttggcgtaatcatggtcatagctgtttcctgtgtgaaattgttatccgctcacaattccacacaacatacgagccggaagcataaagtgtaaagcctggggtgcctaatgagtgagctaactcacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtctgacgcgccctgtagcggcgcattaagcgcggcgggtgtggtggttacgcgcagcgtgaccgctacacttgccagcgccctagcgcccgctcctttcgctttcttcccttcctttctcgccacgttcgccggctttccccgtcaagctctaaatcgggggctccctttagggttccgatttagtgctttacggcacctcgaccccaaaaaacttgattagggtgatggttcacgtagtgggccatcgccctgatagacggtttttcgccctttgacgttggagtccacgttctttaatagtggactcttgttccaaactggaacaacactcaaccctatctcggtctattcttttgatttataagggattttgccgatttcggcctattggttaaaaaaatgagctgatttaacaaaaatttaacgcgaattttaacaaaatattaacgcttacaatttccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcgggcctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcgattaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgacggccagtgagcgcgcctagttattaatagtaatcaattacggggtcattagttcatagcccatatatggagttccgcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatggtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagctctctggctaactagagaacccactgcttactggcttatcgaaattaatacgactcactatagggagacccaagctggctagattaagcgtctgatgagtccgtgaggacgaaacccggcgtaccgggtc'.upper()

    # concatenated to mimic circular plasmid
    ref_seq = 2*(upstream_MCS + 'N'*insert_len + downstream_MCS)
    ref_name = 'SBARRO'
    ref_path = f'{outpath}/.{exp_name}.concat.ref.fa'

    # write ref fasta
    with open(ref_path, 'w') as fh:
        fh.write(f'>SBARRO\n{ref_seq}\n')

    # return ref length for downstream use in output plots
    return len(upstream_MCS + 'N'*insert_len + downstream_MCS)

def process_ref(outpath: str, plasmid_path: str, insert_len: int) -> int:
    # to hold full ref sequence before concatenation
    exp_name = os.path.basename(outpath)
    full_seq = ''
    with open(os.path.normpath(plasmid_path), 'r') as fh:
        with open(f'{outpath}/.{exp_name}.concat.ref.fa', 'w') as out:
            for i,line in enumerate(fh):
                if i == 0:
                    out.write(line)
                else:
                    full_seq+=line
            # remove any whitespace
            full_seq = ''.join(full_seq.split())
            full_seq_concat = 2 * full_seq
            out.write(full_seq_concat.strip() + '\n')
    # return ref length for downstream use in output plots
    return len(full_seq) + insert_len

def _append_doubled_ref(ref_fa_path: str, dest_path: str) -> None:
    """Read a single-record FASTA, concatenate the sequence to itself, and append to dest."""
    with open(ref_fa_path, 'r') as fh:
        header = fh.readline().rstrip()
        seq = ''.join(line.strip() for line in fh)
    with open(dest_path, 'a') as out:
        out.write(f'{header}\n{seq}{seq}\n')

def mm2_align(outpath: str, trimmed_reads_path: str, ap_flag: bool, is_default_plasmid: bool) -> Dict[str, int]:
    exp_name = os.path.basename(outpath)
    plasmid_path = f'{outpath}/.{exp_name}.concat.ref.fa'

    output_bam_file = f'{outpath}/{exp_name}.bam'
    output_sorted_bam_file = f'{outpath}/{exp_name}.sorted.bam'
    output_fa_aligned = f'{outpath}/{exp_name}.aligned.fa'
    output_fa_unaligned = f'{outpath}/{exp_name}.unaligned.fa.gz'
    rev_complement_cmd = 'while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done'

    # E. coli ref:
    ecoli_ref_fa = files('longbarcodeqc.plasmids').joinpath('Escherichia_coli_gca_001606525.ASM160652v1_.dna_rm.toplevel.fa')
    ecoli_ref_string = 'STEC1686_contig' # to grab ecoli reads from bam (all ecoli seqs start with this string)
    ecoli_fa_aligned = f'{outpath}/{exp_name}.ecoli.aligned.fa'
    # Write tmp combined reference of plasmid + ecoli (for quantifying ecoli reads)
    # (without using combined reference, many reads will "align" to both plasmid and ecoli)
    os.system(f'cat {ecoli_ref_fa} {plasmid_path} >{outpath}/.tmp.ref.fa')

    # AP reference strings (match FASTA headers)
    ap_kan_ref_string = 'AP-Kan'
    ap_amp_ref_string = 'AP-Amp'

    # Determine which AP references to add
    # Default plasmid: AP-Amp is already the target, only add AP-Kan
    # -a flag: custom plasmid is target, add both AP-Amp and AP-Kan
    add_ap_kan = is_default_plasmid or ap_flag
    add_ap_amp = ap_flag  # only when -a flag with custom/SBARRO plasmid

    if add_ap_kan:
        ap_kan_ref_fa = files('longbarcodeqc.plasmids').joinpath('AP-Kan.fa')
        output_ap_kan_aligned = f'{outpath}/{exp_name}.AP-Kan.aligned.fa'
        _append_doubled_ref(ap_kan_ref_fa, f'{outpath}/.tmp.ref.fa')

    if add_ap_amp:
        ap_amp_ref_fa = files('longbarcodeqc.plasmids').joinpath('AP-Amp.fa')
        output_ap_amp_aligned = f'{outpath}/{exp_name}.AP-Amp.aligned.fa'
        _append_doubled_ref(ap_amp_ref_fa, f'{outpath}/.tmp.ref.fa')

    # Extract plasmid reference name from fasta, for pulling out plasmid reads from the combined alignment
    with open(plasmid_path, 'r') as fh:
        ref_name = str(fh.readline().strip().split(' ')[0][1:])
    if not os.path.exists(output_fa_unaligned):
        # align reads with mm2, convert to bam output
        print(f'Aligning {exp_name} using minimap2...')
        os.system(f'minimap2 -ax map-ont {outpath}/.tmp.ref.fa {trimmed_reads_path} \
                  --secondary=no -t 3 2>>{outpath}/Log.txt | samtools view -b -h -@ 2 >{output_bam_file}')
        # sort and index bam
        os.system(f'samtools sort {output_bam_file} >{output_sorted_bam_file} 2>>{outpath}/Log.txt')
        os.system(f'samtools index {output_sorted_bam_file}')
        # rm tmp files
        os.system(f'rm {output_bam_file} && rm {outpath}/.tmp.ref.fa && rm {plasmid_path}')

        print(f'Alignment finished. Flipping reverse reads and writing fasta outputs')

        # write positive strand alignments to fasta (only reads aligned to reference name)
        # SAM flags used:
        ## 2048: supplementary alignments
        ## 4:    unmapped
        ## 16:   reverse strand
        ## -q 2: exclude MAPQ <= 1 (ambiguous alignments)
        os.system(f'samtools view -F 2048 -F 4 -F 16 -q 2 -h {output_sorted_bam_file} \'{ref_name}\' | \
                  samtools fasta - >{output_fa_aligned} 2>>{outpath}/Log.txt')
        # reverse complement negative strand alignments and append to same file
        os.system(f'samtools view -F 2048 -f 16 -q 2 -h {output_sorted_bam_file} \'{ref_name}\' | \
                  samtools fasta - 2>>{outpath}/Log.txt | {rev_complement_cmd} >>{output_fa_aligned}')
        os.system(f'gzip {output_fa_aligned}') # compress into .fa.gz

        # write AP-Kan reads to file
        if add_ap_kan:
            os.system(f'samtools view -F 2048 -F 4 -F 16 -q 2 -h {output_sorted_bam_file} \'{ap_kan_ref_string}\' | \
                      samtools fasta - >{output_ap_kan_aligned} 2>>{outpath}/Log.txt')
            os.system(f'samtools view -F 2048 -f 16 -q 2 -h {output_sorted_bam_file} \'{ap_kan_ref_string}\' | \
                      samtools fasta - 2>>{outpath}/Log.txt | {rev_complement_cmd} >>{output_ap_kan_aligned}')
            os.system(f'gzip {output_ap_kan_aligned}')

        # write AP-Amp reads to file (only when -a flag with custom/SBARRO plasmid)
        if add_ap_amp:
            os.system(f'samtools view -F 2048 -F 4 -F 16 -q 2 -h {output_sorted_bam_file} \'{ap_amp_ref_string}\' | \
                      samtools fasta - >{output_ap_amp_aligned} 2>>{outpath}/Log.txt')
            os.system(f'samtools view -F 2048 -f 16 -q 2 -h {output_sorted_bam_file} \'{ap_amp_ref_string}\' | \
                      samtools fasta - 2>>{outpath}/Log.txt | {rev_complement_cmd} >>{output_ap_amp_aligned}')
            os.system(f'gzip {output_ap_amp_aligned}')

        # write unaligned reads to file
        os.system(f'samtools view -f 4 -h {output_sorted_bam_file} | samtools fasta - 2>>{outpath}/Log.txt | \
                  gzip >{output_fa_unaligned}')

        # write e coli reads to file
        os.system(f'samtools view -H {output_sorted_bam_file} >{outpath}/.tmp.ecoli.sam')
        os.system(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                    awk \'$3 ~ /{ecoli_ref_string}/\' >>{outpath}/.tmp.ecoli.sam')
        os.system(f'samtools fasta {outpath}/.tmp.ecoli.sam >{ecoli_fa_aligned} 2>>{outpath}/Log.txt')
        os.system(f'gzip {ecoli_fa_aligned} && rm {outpath}/.tmp.ecoli.sam')

        # store counts for aligned reads, unaligned reads, ecoli reads
        summary_dict = {}
        # Label target as 'AP-Amp' when using default plasmid, otherwise 'Target plasmid'
        target_label = 'AP-Amp' if is_default_plasmid else 'Target plasmid'
        summary_dict[target_label] = int(os.popen(f'samtools view -F 4 -F 2048 -q 2 {output_sorted_bam_file} \'{ref_name}\' | \
                                                    wc -l').read())
        if add_ap_kan:
            summary_dict['AP-Kan'] = int(os.popen(f'samtools view -F 4 -F 2048 -q 2 {output_sorted_bam_file} \'{ap_kan_ref_string}\' | \
                                                        wc -l').read())
        if add_ap_amp:
            summary_dict['AP-Amp'] = int(os.popen(f'samtools view -F 4 -F 2048 -q 2 {output_sorted_bam_file} \'{ap_amp_ref_string}\' | \
                                                        wc -l').read())
        summary_dict['Ambiguous alignments'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                                                    awk \'$5 <= 1\' | wc -l').read())
        summary_dict['Unaligned'] = int(os.popen(f'samtools view -f 4 {output_sorted_bam_file} | wc -l').read())
        summary_dict['E. coli'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                                                    awk \'$3 ~ /{ecoli_ref_string}/\' | wc -l').read())

        read_total = sum(summary_dict.values())
        print(f'\nAlignment stats ({read_total} total reads):')
        print(f'{summary_dict[target_label]} {target_label} reads aligned')
        if add_ap_kan: print(f'{summary_dict["AP-Kan"]} AP-Kan reads aligned')
        if add_ap_amp: print(f'{summary_dict["AP-Amp"]} AP-Amp reads aligned')
        print(f'{summary_dict["E. coli"]} bacterial reads aligned')
        print(f'{summary_dict["Ambiguous alignments"]} ambiguous alignments (MAPQ <= 1)')
        print(f'{summary_dict["Unaligned"]} unaligned reads\n')
        return summary_dict
