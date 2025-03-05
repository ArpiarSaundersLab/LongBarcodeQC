import parasail
import pandas as pd
import numpy as np
import os
from importlib.resources import files
from tqdm import tqdm

def barcode_scores(outpath: str, barcode_path: str, flanks_path, 
                   insert_len: int, enzymes, a31, SBARRO) -> None:
    exp_name = os.path.basename(outpath)
    # loading reads from the experiment
    long_reads = parasail.sequences_from_file(f'{outpath}/{exp_name}.aligned.fa.gz')
    
    # loading possible barcode sequences
    barcodes = parasail.sequences_from_file(barcode_path)
    bc_name_seq = [(barcodes[i].name.decode(), barcodes[i].seq.decode()) for i in range(len(barcodes))]

    # define restriction sites
    restriction_sites = [
        ('RE_AvrII','CCTAGG'),
        ('RE_KpnI','GGTACC'),
        ('RE_PciI','ACATGT'),
        ('RE_SpeI','ACTAGT'),
        ('RE_PlutI','GGCGCC'),
        ('RE_NotI','GCGGCCGC'),
        ('RE_TspmI','CCCGGG'),
        ('RE_MreI','CGCCGGCG'),
        ('RE_MauBI','CGCGCGCG'),
        ('RE_BsiWI','CGTACG')
    ]

    if enzymes: # add user defined enzymes if provided
        with open(os.path.normpath(enzymes), 'r') as fh:
            for rs in fh:
                name = str(rs.strip().split(',')[0])
                seq = str(rs.strip().split(',')[1])
                restriction_sites.append((f'RE_{name}',seq))

    # build common MCS flanking region to anchor to the barcode region for a more precise alignment
    # These flanking regions are common for all c.18 derivatives. Add NNN sequence of user specified
    # length between them.
    if SBARRO:
        SBARRO_flanks_path = str(files('barcodezy.plasmids').joinpath('SBARRO_flanks.fa'))
        flanks = parasail.sequences_from_file(SBARRO_flanks_path)
        MCS_left_flank = str(flanks[0].seq.decode())
        MCS_right_flank = str(flanks[1].seq.decode())
        # insert N sequence of expected insert length between flanks
        MCS_flank = MCS_left_flank + 'N'*insert_len + MCS_right_flank  
    else:
        f_path = str(os.path.normpath(flanks_path))
        flanks = parasail.sequences_from_file(f_path)
        if len(flanks) != 2:
            raise ValueError('Error: Flank file must contain 2 sequences: '
                             'Upstream flank followed by downstream flank. '
                             'Recommended that each flank be ~75 bp.')
        MCS_left_flank = str(flanks[0].seq.decode())
        MCS_right_flank = str(flanks[1].seq.decode())
        # insert N sequence of expected insert length between flanks
        MCS_flank = MCS_left_flank + 'N'*insert_len + MCS_right_flank  

    # main dict to hold lists (columns) of read statistics
    df_dict = {
        'seq_id': [],
        'seq_len': [],
        'read_type': [],
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
    MCS_query = parasail.profile_create_16(MCS_flank, user_matrix)

    # iterate through reads and extract output stats
    print('Processing target plasmid reads...')
    for i,longread in enumerate(tqdm(long_reads, desc='Aligning barcodes (target plasmid)',
                                     bar_format="{desc}: |{bar}| {percentage:3.0f}% ({n} reads)",
                                     leave=False)):
        # avoid weird sequence objects possibly causing segmentation fault
        if i > (len(long_reads) - 1):
            # the final sequence object is empty for some reason
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
        if (MCS_len < len(MCS_flank)):
            df_dict['MCS_score'].append(MCS_score)
            df_dict['read_type'].append('plasmid_failed_anchor')
            df_dict['seq_len'].append(read_len)
            df_dict['MCS_len'].append(MCS_len)
            df_dict['seq_id'].append(read_id)
            [df_dict[bc[0]].append(np.nan) for bc in bc_name_seq]
            [df_dict[rs[0]].append(rs[1] in longread.seq.decode()) for rs in restriction_sites]
            continue

        # get barcode scores
        bc_query = parasail.profile_create_16(longread_doubled[MCS_start:MCS_end], user_matrix)
        for bc in bc_name_seq:
            bc_align = parasail.sg_qx_striped_profile_16(bc_query, bc[1], 5, 1)
            df_dict[bc[0]].append(bc_align.score)

        # append df values
        df_dict['MCS_score'].append(MCS_score)
        df_dict['MCS_len'].append(MCS_len)
        # get restriction site booleans
        [df_dict[rs[0]].append(rs[1] in longread_doubled[MCS_start:MCS_end]) for rs in restriction_sites]
        df_dict['read_type'].append('plasmid')
        df_dict['seq_len'].append(read_len)
        df_dict['seq_id'].append(read_id)
    tqdm.write('Done\n')

    #### 
    # call a31 barcodes if user flags a31
    ####
    if a31:
        print('Processing a.31 reads...')
        long_reads_a31 = parasail.sequences_from_file(f'{outpath}/{exp_name}.a31.aligned.fa.gz')

        a31_flanks_path = str(files('barcodezy.plasmids').joinpath('a.31_flanks.fa'))
        flanks = parasail.sequences_from_file(a31_flanks_path)
        MCS_left_flank = str(flanks[0].seq.decode())
        MCS_right_flank = str(flanks[1].seq.decode())
        # insert N sequence of expected insert length between flanks
        MCS_flank = MCS_left_flank + 'N'*insert_len + MCS_right_flank
            
        #custom align matrix
        user_matrix = parasail.matrix_create("ACGT", match=5, mismatch=-3)
        MCS_query = parasail.profile_create_16(MCS_flank, user_matrix)

        for k,longread in enumerate(tqdm(long_reads_a31, desc='Aligning barcodes (a31)',
                                         bar_format="{desc}: |{bar}| {percentage:3.0f}% ({n} reads)",
                                         leave=False)):
            # avoid weird sequence objects possibly causing segmentation fault
            if k > (len(long_reads_a31) - 1):
                # again, the final sequence object is empty for some reason
                continue

            read_id = longread.name.decode()
            read_len = len(longread)

            # double the sequence to mimic linear sequence
            longread_doubled = longread.seq.decode() * 2

            align_MCS = parasail.sg_dx_striped_profile_16(MCS_query, longread_doubled, 5, 1)

            MCS_start = align_MCS.end_ref - len(MCS_flank)
            MCS_end = align_MCS.end_ref
            MCS_score = align_MCS.score
            MCS_len = len(longread_doubled[MCS_start:MCS_end])
        
            # catch reads where MCS flank isn't fully covered
            if (MCS_len < len(MCS_flank)):
                df_dict['MCS_score'].append(MCS_score)
                df_dict['read_type'].append('a31_failed_anchor')
                df_dict['seq_len'].append(read_len)
                df_dict['MCS_len'].append(MCS_len)
                df_dict['seq_id'].append(read_id)
                [df_dict[bc[0]].append(np.nan) for bc in bc_name_seq]
                [df_dict[rs[0]].append(rs[1] in longread.seq.decode()) for rs in restriction_sites]
                continue

            # get barcode scores
            bc_query = parasail.profile_create_16(longread_doubled[MCS_start:MCS_end], user_matrix)
            for bc in bc_name_seq:
                bc_align = parasail.sg_qx_striped_profile_16(bc_query, bc[1], 5, 1)
                df_dict[bc[0]].append(bc_align.score)

            # append df values
            df_dict['MCS_score'].append(MCS_score)
            df_dict['read_type'].append('a31')
            df_dict['seq_len'].append(read_len)
            df_dict['MCS_len'].append(MCS_len)
            df_dict['seq_id'].append(read_id)
            # get restriction site booleans
            [df_dict[rs[0]].append(rs[1] in longread_doubled[MCS_start:MCS_end]) for rs in restriction_sites]
    tqdm.write('Done\n')
    
    ### ecoli reads
    print('Processing E. coli reads...')
    long_reads_ecoli = parasail.sequences_from_file(f'{outpath}/{exp_name}.ecoli.aligned.fa.gz') 
    for z,longread in enumerate(long_reads_ecoli):
        # avoid weird sequence objects possibly causing segmentation fault
        if z > (len(long_reads_ecoli) - 1):
            # again, the final sequence object is empty for some reason
            continue

        read_id = longread.name.decode()
        read_len = len(longread)
        for bc in bc_name_seq:
            df_dict[bc[0]].append(np.nan)

        # append df values
        df_dict['MCS_score'].append(np.nan)
        df_dict['read_type'].append('ecoli')
        df_dict['seq_len'].append(read_len)
        df_dict['MCS_len'].append(np.nan)
        df_dict['seq_id'].append(read_id)
        # get restriction site booleans
        [df_dict[rs[0]].append(np.nan) for rs in restriction_sites]
    print('Done\n')

    print(f'Finished read processing and barcode alignments\nGenerating outputs...')

    output_df = pd.DataFrame.from_dict(df_dict)
    return output_df
