import sys, os
from importlib.resources import files
from typing import Dict

def rename_reads(raw_read_file, exp_name) -> None:
    print(f'Writing file...')
    with open(raw_read_file, 'r') as reads:
        tmp_rewrite_file = f'{raw_read_file[0:-5]}_renamed.fastq'
        with open(tmp_rewrite_file, 'w') as out:
            read_count=0
            for i,line in enumerate(reads):
                if i % 4 == 0:
                    read_count+=1
                    id = line.strip().split(' ')
                    out.write(f'@{exp_name}_{read_count} {' '.join(id[1:])}\n')
                if i % 4 == 1:
                    seq = line.strip()
                    out.write(f'{seq}\n')
                if i % 4 == 2:
                    strand = line.strip()
                    out.write(f'{strand}\n')
                if i % 4 == 3:
                    qual = line.strip()
                    out.write(f'{qual}\n')

    os.replace(tmp_rewrite_file, raw_read_file)

def generate_ref(outpath, insert_len):
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

def process_ref(outpath, plasmid_path):
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
            out.write(full_seq_concat.strip())
    # return ref length for downstream use in output plots
    return len(full_seq)

def mm2_align(outpath, trimmed_reads_path, a31) -> Dict[str, int]:
    exp_name = os.path.basename(outpath)
    plasmid_path = f'{outpath}/.{exp_name}.concat.ref.fa'

    output_bam_file = f'{outpath}/{exp_name}.bam'
    output_sorted_bam_file = f'{outpath}/{exp_name}.sorted.bam'
    output_fa_aligned = f'{outpath}/{exp_name}.aligned.fa'
    output_fa_unaligned = f'{outpath}/{exp_name}.unaligned.fa.gz'
    rev_complement_cmd = 'while read L; do  echo $L; read L; echo "$L" | rev | tr "ATGC" "TACG" ; done'

    # E. coli ref:
    ecoli_ref_fa = files('barcodezy.plasmids').joinpath('Escherichia_coli_gca_001606525.ASM160652v1_.dna_rm.toplevel.fa')
    ecoli_ref_string = 'STEC1686_contig' # to grab ecoli reads from bam (all ecoli seqs start with this string)
    ecoli_fa_aligned = f'{outpath}/{exp_name}.ecoli.aligned.fa'
    # Write tmp combined reference of plasmid + ecoli (for quantifying ecoli reads)
    # (without using combined reference, many reads will "align" to both plasmid and ecoli)
    os.system(f'cat {ecoli_ref_fa} {plasmid_path} >{outpath}/.tmp.ref.fa')

    # if a31 flag turned on, add a31 to reference
    if a31:
        # a.31 ref:
        output_a31_aligned = f'{outpath}/{exp_name}.a31.aligned.fa'
        a31_ref_fa = files('barcodezy.plasmids').joinpath('a.31.fa')
        a31_ref_string = 'a.31_pUC57-Mini_STICR_6_v2_circ'
        os.system(f'cat {a31_ref_fa} >>{outpath}/.tmp.ref.fa')

    # Extract plasmid reference name from fasta, for pulling out plasmid reads from the combined alignment
    with open(plasmid_path, 'r') as fh:
        ref_name = str(fh.readline().strip().split(' ')[0][1:])
    if not os.path.exists(output_fa_unaligned):
        # align reads with mm2, convert to bam output
        print(f'\nAligning {exp_name} with minimap2')
        os.system(f'minimap2 -ax map-ont {outpath}/.tmp.ref.fa {trimmed_reads_path} \
        --secondary=no -t 3 | samtools view -b -h -@ 2 >{output_bam_file}')
        # sort and index bam
        os.system(f'samtools sort {output_bam_file} >{output_sorted_bam_file}')
        os.system(f'samtools index {output_sorted_bam_file}')
        # rm tmp files
        os.system(f'rm {output_bam_file} && rm {outpath}/.tmp.ref.fa && rm {plasmid_path}')

        print(f'\nAlignment finished. Flipping reverse reads and writing fasta outputs')

        # write positive strand alignments to fasta (only reads aligned to reference name)
        # SAM flags used:
        ## 2048: supplementary alignments
        ## 4:    unmapped
        ## 16:   reverse strand
        os.system(f'samtools view -F 2048 -F 4 -F 16 -h {output_sorted_bam_file} \'{ref_name}\' | \
                  samtools fasta - > {output_fa_aligned}')
        # reverse complement negative strand alignments and append to same file
        os.system(f'samtools view -F 2048 -f 16 -h {output_sorted_bam_file} \'{ref_name}\' | samtools fasta - | \
                  {rev_complement_cmd} >> {output_fa_aligned}')
        os.system(f'gzip {output_fa_aligned}') # compress into .fa.gz

        # write a.31 reads to file if a31 flag on:
        if a31:
            os.system(f'samtools view -F 2048 -F 4 -F 16 -h {output_sorted_bam_file} \'{a31_ref_string}\' | \
                      samtools fasta - > {output_a31_aligned}')
            # reverse complement negative strand alignments and append to same file
            os.system(f'samtools view -F 2048 -f 16 -h {output_sorted_bam_file} \'{a31_ref_string}\' | samtools fasta - | \
                      {rev_complement_cmd} >> {output_a31_aligned}')
            os.system(f'gzip {output_a31_aligned}') # compress into .fa.gz

        # write unaligned reads to file
        os.system(f'samtools view -f 4 -h {output_sorted_bam_file} | samtools fasta - | gzip > {output_fa_unaligned}')

        # write e coli reads to file
        os.system(f'samtools view -H {output_sorted_bam_file} >{outpath}/.tmp.ecoli.sam')
        os.system(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                    awk \'$3 ~ /{ecoli_ref_string}/\' >>{outpath}/.tmp.ecoli.sam')
        os.system(f'samtools fasta {outpath}/.tmp.ecoli.sam >{ecoli_fa_aligned}')
        os.system(f'gzip {ecoli_fa_aligned} && rm {outpath}/.tmp.ecoli.sam')

        
        # store counts for aligned reads, unaligned reads, ecoli reads
        summary_dict = {}
        summary_dict['Target plasmid'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} \'{ref_name}\' | \
                                                    wc -l').read())
        if a31:
            summary_dict['a.31'] = int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} \'{a31_ref_string}\' | \
                                                        wc -l').read())
        summary_dict['Unaligned'] = int(os.popen(f'samtools view -f 4 {output_sorted_bam_file} | wc -l').read())
        summary_dict['E. coli']= int(os.popen(f'samtools view -F 4 -F 2048 {output_sorted_bam_file} | \
                                                    awk \'$3 ~ /{ecoli_ref_string}/\' | wc -l').read()) 

        read_total = sum(summary_dict.values())
        print(f'\nAlignment stats ({read_total} total reads):')
        print(f'{summary_dict["Target plasmid"]} plasmid reads aligned')
        if a31: print(f'{summary_dict["a.31"]} a.31 reads aligned')
        print(f'{summary_dict["E. coli"]} bacterial reads aligned')
        print(f'{summary_dict["Unaligned"]} undetermined reads\n')
        return summary_dict