import parasail
import pandas as pd

def barcode_scores(reads_file: str, barcode_path: str, insert_len: int) -> None:
    # loading reads from the experiment
    long_reads = parasail.sequences_from_file(reads_file)
    
    # loading possible barcode sequences
    barcodes = parasail.sequences_from_file(barcode_path)
    bc_name_seq = [(barcodes[i].name.decode(), barcodes[i].seq.decode()) for i in range(len(barcodes))]

    # define restriction sites
    restriction_sites = [
        ('AvrII_site1','CCTAGG'),
        ('KpnI_site2','GGTACC'),
        ('Pcil_site3','ACATGT'),
        ('SpeI_site4_1','ACTAGT'),
        ('PlutI_site4_2','GGCGCC'),
        ('NotI_site4_3','GCGGCCGC'),
        ('TspmI_c28','CCCGGG'),
        ('MreI','CGCCGGCG'),
        ('MauBI','CGCGCGCG'),
        ('BsiWI','CGTACG')
    ]

    # build common MCS flanking region to anchor to the barcode region for a more precise alignment
    MCS_left_flank = ('CAAAAAAGAAGAGAAAGGTAGATCCAAAAAAGAAGAGAAAGGTAGATCCAAAAAAGAAGAGAAAGGTAGGATCCA'
                      'GATAACTGATCATAATCAGCCATACCACATTTGGAATTCGTCGAGGGGGCCGCCGGCGAGGCCGGCCATGAATTC')
    MCS_right_flank = ('ACCGGTTACTTGCGATCGCAATCGACGCGTTCCTGCAGGACGCGCGCGGCTAGACATGAAAAAAACTAACACTCC'
                       'TCCGGTACCGCCACCATGGACGTGGTGAATCAGCTGGTGGCTGGGGGTCAGTTCCGGGTGGTCAAGGAGCCCCTT')
    # insert N sequence of expected insert length between flanks
    MCS_flank = MCS_left_flank + 'N'*insert_len + MCS_right_flank  

    # main dict to hold lists (columns) of read statistics
    df_dict = {
        'seq_id': [],
        'seq_len': [],
        'plasmid_aligned': [],
        'MCS_score': [],
        'MCS_len': []
    }

    # initialize output lists for each restriction site in df_dict
    for rs in restriction_sites:
        df_dict[rs[0]] = []
    # initialize output lists for each barcode in df_dict
    for bc in bc_name_seq:
        df_dict[bc[0]] = []

    #custom align matrix
    user_matrix = parasail.matrix_create("ACGT", match=5, mismatch=-3)
    # define alignment profile of MCS_flank before loop (should improve speed)
    MCS_query = parasail.profile_create_16(MCS_flank, user_matrix)
    # iterate through reads and extract output stats
    j=1000
    for i,longread in enumerate(long_reads):
        # progress counter per 1000 reads
        if i % 1000 == 0 and i > 0:
            print(f'{j} reads processed')
            j+=1000

        # avoid weird sequence objects possibly causing segmentation fault
        if i > (len(long_reads) - 2):
            #print('Warning: skipped empty parasail seq object')
            continue

        read_id = longread.name.decode()
        read_len = len(longread)

        # anchoring alignment; double the sequence to mimic linear sequence
        longread_doubled = longread.seq.decode() * 2
        align_MCS = parasail.sg_dx_striped_profile_16(MCS_query, longread_doubled, 5, 1)

        MCS_start = align_MCS.end_ref - len(MCS_flank)
        MCS_end = align_MCS.end_ref
        MCS_score = align_MCS.score
        MCS_len = len(longread_doubled[MCS_start:MCS_end])
        
        # catch reads where MCS flank isn't fully covered
        if MCS_len < len(MCS_flank):
            df_dict['MCS_score'].append(0)
            df_dict['plasmid_aligned'].append(False)
            df_dict['seq_len'].append(read_len)
            df_dict['MCS_len'].append(0)
            df_dict['seq_id'].append(read_id)
            [df_dict[bc[0]].append(0) for bc in bc_name_seq]
            [df_dict[rs[0]].append(rs[1] in longread.seq.decode()) for rs in restriction_sites]
            continue

        # get barcode scores
        bc_query = parasail.profile_create_16(longread_doubled[MCS_start:MCS_end], user_matrix)
        for bc in bc_name_seq:
            bc_align = parasail.sg_qx_striped_profile_16(bc_query, bc[1], 5, 1)
            df_dict[bc[0]].append(bc_align.score)

        # append df values
        df_dict['MCS_score'].append(MCS_score)
        df_dict['plasmid_aligned'].append(True)
        df_dict['seq_len'].append(read_len)
        df_dict['MCS_len'].append(MCS_len)
        df_dict['seq_id'].append(read_id)
        # get restriction site booleans
        [df_dict[rs[0]].append(rs[1] in longread_doubled[MCS_start:MCS_end]) for rs in restriction_sites]

    print(f'Finished read processing\n Writing csv output\n Processed {i} reads')

    output_df = pd.DataFrame.from_dict(df_dict)
    return output_df