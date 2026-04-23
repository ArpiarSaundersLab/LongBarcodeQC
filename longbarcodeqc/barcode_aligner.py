from typing import Optional

import os

import numpy as np
import pandas as pd
import parasail
from importlib.resources import files
from tqdm import tqdm

_GAP_OPEN = 5
_GAP_EXTEND = 2


def _load_flanks(flanks_path: str, insert_len: int) -> tuple[str, str, int]:
    """Load MCS flanking sequences; return (left_flank, right_flank, full_flank_len)."""
    flanks = parasail.sequences_from_file(flanks_path)
    if len(flanks) != 2:
        raise ValueError(
            'Flank file must contain 2 sequences: upstream flank followed by downstream flank. '
            'Recommended that each flank be ~75 bp.'
        )
    left = flanks[0].seq.decode()
    right = flanks[1].seq.decode()
    return left, right, len(left) + insert_len + len(right)


def _find_mcs(seq: str, query_left, query_right, right_flank_len: int) -> tuple[int, int, int, int]:
    """Align MCS flanks to a sequence; return (start, end, left_score, right_score)."""
    left_result = parasail.sg_striped_profile_16(query_left, seq, _GAP_OPEN, _GAP_EXTEND)
    right_result = parasail.sg_striped_profile_16(query_right, seq, _GAP_OPEN, _GAP_EXTEND)
    return (
        left_result.end_ref,
        right_result.end_ref - right_flank_len,
        left_result.score,
        right_result.score,
    )


def _process_read(
    longread,
    read_type: str,
    bc_name_seq: list[tuple[str, str]],
    restriction_sites: list[tuple[str, str]],
    query_left,
    query_right,
    mcs_flank_len: int,
    right_flank_len: int,
    user_matrix,
) -> dict:
    """Process a single read; return a row dict of stats."""
    read_id = longread.name.decode()
    read_len = len(longread)
    seq = longread.seq.decode()
    doubled = seq * 2

    start, end, left_score, right_score = _find_mcs(doubled, query_left, query_right, right_flank_len)

    # fallback: try the middle window of the doubled sequence
    if start >= end or (end - start) > 2 * mcs_flank_len:
        doubled = doubled[read_len // 2: (read_len // 2) + read_len]
        start, end, left_score, right_score = _find_mcs(doubled, query_left, query_right, right_flank_len)

    failed = start >= end or (end - start) > 2 * mcs_flank_len
    mcs_seq = doubled[start:end]

    row = {
        'seq_id': read_id,
        'seq_len': read_len,
        'read_type': f'{read_type}_failed_anchor' if failed else read_type,
        'left_anchor_score': left_score,
        'right_anchor_score': right_score,
        'MCS_len': len(mcs_seq),
        'MCS_seq': mcs_seq,
    }

    if failed:
        for bc_name, _ in bc_name_seq:
            row[bc_name] = np.nan
        # check restriction sites on full read (MCS region not reliably extracted)
        for rs_name, rs_seq in restriction_sites:
            row[rs_name] = rs_seq in seq
    else:
        bc_query = parasail.profile_create_16(mcs_seq, user_matrix)
        for bc_name, bc_seq in bc_name_seq:
            row[bc_name] = parasail.sg_striped_profile_16(bc_query, bc_seq, _GAP_OPEN, _GAP_EXTEND).score
        for rs_name, rs_seq in restriction_sites:
            row[rs_name] = rs_seq in mcs_seq

    return row


def _process_batch(
    reads,
    read_type: str,
    bc_name_seq: list[tuple[str, str]],
    restriction_sites: list[tuple[str, str]],
    left_flank: str,
    right_flank: str,
    mcs_flank_len: int,
    user_matrix,
    desc: str,
) -> list[dict]:
    """Process a batch of reads; return a list of row dicts."""
    query_left = parasail.profile_create_16(left_flank, user_matrix)
    query_right = parasail.profile_create_16(right_flank, user_matrix)
    right_flank_len = len(right_flank)

    rows = []
    for i, longread in enumerate(tqdm(reads, desc=desc, bar_format="{desc}: |{bar}| {percentage:3.0f}% ({n} reads)", leave=False)):
        # parasail sequences iterator yields one extra invalid object beyond len()
        if i > (len(reads) - 1):
            continue
        rows.append(_process_read(
            longread, read_type, bc_name_seq, restriction_sites,
            query_left, query_right, mcs_flank_len, right_flank_len, user_matrix,
        ))
    tqdm.write('Done\n')
    return rows


def barcode_scores(
    outpath: str,
    barcode_path: str,
    flanks_path: Optional[str],
    insert_len: int,
    enzymes: Optional[str],
    ap_flag: bool,
    is_default_plasmid: bool,
    SBARRO: bool,
) -> pd.DataFrame:
    """Compute barcode alignment scores and read statistics.

    Returns a DataFrame with one row per read and columns for scores/metadata.
    """
    exp_name = os.path.basename(outpath)

    barcodes = parasail.sequences_from_file(barcode_path)
    bc_name_seq = [(barcodes[i].name.decode(), barcodes[i].seq.decode()) for i in range(len(barcodes))]

    restriction_sites = [
        ('RE_AvrII', 'CCTAGG'),
        ('RE_KpnI',  'GGTACC'),
        ('RE_PciI',  'ACATGT'),
        ('RE_SpeI',  'ACTAGT'),
        ('RE_NotI',  'GCGGCCGC'),
        ('RE_AsiSI',  'GCGATCGC'),
        ('RE_BsiWI', 'CGTACG'),
        ('RE_MreI',  'CGCCGGCG'),
        ('RE_FseI', 'GGCCGGCC'),
        ('RE_Sbf1', 'CCTGCAGG'),
        ('RE_MauBI', 'CGCGCGCG')
    ]

    if enzymes:
        with open(os.path.normpath(enzymes)) as fh:
            for rs in fh:
                name, seq = rs.strip().split(',')[:2]
                restriction_sites.append((f'RE_{name}', seq))

    if SBARRO:
        flanks_path = str(files('longbarcodeqc.plasmids').joinpath('SBARRO_flanks.fa'))
    left_flank, right_flank, mcs_flank_len = _load_flanks(flanks_path, insert_len)

    user_matrix = parasail.matrix_create("ACGT", match=2, mismatch=-1)

    # process target plasmid reads
    target_read_type = 'AP-Amp' if is_default_plasmid else 'plasmid'
    print(f'Processing {target_read_type} reads...')
    long_reads = parasail.sequences_from_file(f'{outpath}/{exp_name}.aligned.fa.gz')
    rows = _process_batch(
        long_reads, target_read_type, bc_name_seq, restriction_sites,
        left_flank, right_flank, mcs_flank_len, user_matrix,
        desc=f'Aligning barcodes ({target_read_type})',
    )

    # process AP reads if applicable
    add_ap_kan = is_default_plasmid or ap_flag
    add_ap_amp = ap_flag

    if add_ap_kan or add_ap_amp:
        ap_flanks_path = str(files('longbarcodeqc.plasmids').joinpath('AP_flanks.fa'))
        ap_left, ap_right, ap_mcs_len = _load_flanks(ap_flanks_path, insert_len)

    if add_ap_kan:
        print('Processing AP-Kan reads...')
        long_reads_ap_kan = parasail.sequences_from_file(f'{outpath}/{exp_name}.AP-Kan.aligned.fa.gz')
        rows += _process_batch(
            long_reads_ap_kan, 'AP-Kan', bc_name_seq, restriction_sites,
            ap_left, ap_right, ap_mcs_len, user_matrix,
            desc='Aligning barcodes (AP-Kan)',
        )

    if add_ap_amp:
        print('Processing AP-Amp reads...')
        long_reads_ap_amp = parasail.sequences_from_file(f'{outpath}/{exp_name}.AP-Amp.aligned.fa.gz')
        rows += _process_batch(
            long_reads_ap_amp, 'AP-Amp', bc_name_seq, restriction_sites,
            ap_left, ap_right, ap_mcs_len, user_matrix,
            desc='Aligning barcodes (AP-Amp)',
        )

    # process e. coli reads (no barcode alignment needed)
    print('Processing E. coli reads...')
    long_reads_ecoli = parasail.sequences_from_file(f'{outpath}/{exp_name}.ecoli.aligned.fa.gz')
    for z, longread in enumerate(long_reads_ecoli):
        if z > (len(long_reads_ecoli) - 1):
            continue
        row: dict = {
            'seq_id': longread.name.decode(),
            'seq_len': len(longread),
            'read_type': 'ecoli',
            'left_anchor_score': np.nan,
            'right_anchor_score': np.nan,
            'MCS_len': np.nan,
            'MCS_seq': np.nan,
        }
        for bc_name, _ in bc_name_seq:
            row[bc_name] = np.nan
        for rs_name, _ in restriction_sites:
            row[rs_name] = np.nan
        rows.append(row)
    print('Done\n')

    print('Finished read processing and barcode alignments\nGenerating outputs...')
    return pd.DataFrame(rows)
