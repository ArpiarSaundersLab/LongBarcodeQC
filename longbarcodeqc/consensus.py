import os
import re
import subprocess
from typing import Optional

import numpy as np
import parasail

from longbarcodeqc.barcode_aligner import _find_mcs, _load_flanks

# only full-length backbones are used: those within this fraction of the most common
# backbone length. Partial reads (common for large plasmids), multimers and plasmids
# with large deletions are left out.
_LENGTH_TOLERANCE = 0.05
# window used to locate the most common backbone length (narrower than the tolerance
# above, or every length within the full-length cluster would tie)
_PEAK_WINDOW = 0.01
_MIN_READS = 5
# each round realigns the reads to the previous consensus, which removes the errors
# carried in from the read used as the starting draft; it usually settles in 3-4 rounds
_MAX_POLISH_ROUNDS = 5
# differences from the reference longer than this are reported as structural and left
# out of the identity, so e.g. one inserted cassette does not mask the small differences
_LARGE_INDEL = 50
# a flank must align to the reference at this fraction of a perfect score for the
# reference backbone to be used for comparison
_REF_FLANK_MIN_SCORE = 0.8
_CS_OP_RE = re.compile(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)')


def _read_fasta(path: str) -> list[tuple[str, str]]:
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                records.append([line[1:].split()[0], []])
            elif line:
                records[-1][1].append(line)
    return [(name, ''.join(chunks).upper()) for name, chunks in records]


def _with_qualities(backbones: list[tuple[str, str]],
                    fastq_path: str) -> list[tuple[str, str, str, int]]:
    """Look up the base qualities of each backbone in the run FASTQ.

    samtools consensus weighs each base by its quality, and calls noticeably fewer
    consensus bases against the read majority with them than without. The backbone is a
    rotation of its read (reverse complemented if the read aligned to the minus strand),
    so it is found as a substring of the oriented read concatenated to itself.

    Returns (name, backbone, qualities, junction distance). The junction is where the read's
    two ends meet, which carries adapter sequence; the distance is how far it sits from the
    nearer backbone end, going either way around the plasmid. Polishing removes adapter
    sequence from the middle of a draft but not from its ends, where reads only soft-clip.
    """
    wanted = dict(backbones)
    complement = str.maketrans('ACGTN', 'TGCAN')
    records = []
    with open(fastq_path) as fh:
        for header in fh:
            seq = fh.readline().strip().upper()
            fh.readline()
            qual = fh.readline().strip()
            name = header[1:].split()[0]
            if name not in wanted:
                continue
            backbone = wanted[name]
            for read, read_qual in ((seq, qual), (seq.translate(complement)[::-1], qual[::-1])):
                i = (read + read).find(backbone)
                if i >= 0:
                    junction = (len(read) - i) % len(read)
                    if junction < len(backbone):
                        distance = min(junction, len(backbone) - junction)
                    else:  # the read was cut inside the MCS
                        distance = min(junction - len(backbone), len(read) - junction)
                    records.append((name, backbone, (read_qual + read_qual)[i:i + len(backbone)],
                                    distance))
                    break
    return records


def _write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n{seq}\n')


def _run(cmd: str, log_path: str) -> None:
    subprocess.run(f'{cmd} 2>>{log_path}', shell=True, check=True)


def _align_sorted(ref_path: str, reads_path: str, bam_path: str, log_path: str) -> None:
    _run(f'minimap2 -ax map-ont --secondary=no -t 3 {ref_path} {reads_path} 2>>{log_path} | '
         f'samtools sort -@ 2 -o {bam_path} -', log_path)
    _run(f'samtools index {bam_path}', log_path)


def _depth(bam_path: str, ref_len: int) -> np.ndarray:
    out = subprocess.run(f'samtools depth -a {bam_path}', shell=True,
                         capture_output=True, text=True, check=True).stdout
    depth = np.zeros(ref_len, dtype=int)
    for line in out.splitlines():
        _, pos, d = line.split('\t')
        depth[int(pos) - 1] = int(d)
    return depth


def _peak_length(lengths: np.ndarray) -> int:
    """Return the backbone length with the most reads within _PEAK_WINDOW of it.

    Full-length reads pile up at the true backbone length, while partial reads spread
    across all shorter lengths, so the densest length marks the full-length plasmid even
    when most reads are partial (where the median would not).
    """
    lengths = np.sort(lengths)
    lo = np.searchsorted(lengths, lengths * (1 - _PEAK_WINDOW), side='left')
    hi = np.searchsorted(lengths, lengths * (1 + _PEAK_WINDOW), side='right')
    return int(lengths[np.argmax(hi - lo)])


def _reference_backbone(ref_seq: str, flanks_path: str, insert_len: int) -> Optional[tuple[str, int]]:
    """Cut the backbone out of the reference plasmid with the same flank anchoring used on reads.

    Returns (backbone, 0-based plasmid position of the backbone's first base), or None if
    the reference does not contain both flanks (e.g. the partial default AP-Amp reference).
    """
    left, right, _ = _load_flanks(flanks_path, insert_len)
    matrix = parasail.matrix_create('ACGT', match=2, mismatch=-1)
    query_left = parasail.profile_create_16(left, matrix)
    query_right = parasail.profile_create_16(right, matrix)
    ref_len = len(ref_seq)
    start, end, left_score, right_score = _find_mcs(ref_seq * 2, query_left, query_right)
    if (start >= end or left_score < _REF_FLANK_MIN_SCORE * 2 * len(left)
            or right_score < _REF_FLANK_MIN_SCORE * 2 * len(right)):
        return None
    return (ref_seq * 4)[end:start + ref_len], end % ref_len


def _compare_to_reference(consensus_path: str, ref_backbone: str, plasmid_offset: int,
                          plasmid_len: int, outpath: str, log_path: str) -> Optional[dict]:
    """Align the consensus to the reference backbone and list the differences."""
    ref_path = f'{outpath}/.tmp.backbone.ref.fa'
    _write_fasta(ref_path, [('reference_backbone', ref_backbone)])
    paf = subprocess.run(f'minimap2 -c --cs -x asm5 {ref_path} {consensus_path} 2>>{log_path}',
                         shell=True, capture_output=True, text=True, check=True).stdout
    os.remove(ref_path)

    hits = [line.split('\t') for line in paf.splitlines()]
    hits = [h for h in hits if 'tp:A:P' in h]
    if not hits:
        return None
    # keep the longest primary alignment
    hit = max(hits, key=lambda h: int(h[10]))
    if hit[4] != '+':
        return None
    q_pos, r_pos = int(hit[2]), int(hit[7])
    cs = next(f[5:] for f in hit[12:] if f.startswith('cs:Z:'))

    variants = []
    matches = 0
    small_diff_bases = 0
    for op in _CS_OP_RE.findall(cs):
        kind, body = op[0], op[1:]
        if kind == ':':
            matches += int(body)
            q_pos += int(body)
            r_pos += int(body)
            continue
        if kind == '*':
            ref_base, cons_base = body[0].upper(), body[1].upper()
            # an N in the consensus is a no-call, not a difference from the reference
            if cons_base != 'N':
                variants.append(('Substitution', r_pos, q_pos, ref_base, cons_base))
                small_diff_bases += 1
            q_pos += 1
            r_pos += 1
        elif kind == '+':
            variants.append(('Insertion', r_pos, q_pos, '', body.upper()))
            q_pos += len(body)
        else:
            variants.append(('Deletion', r_pos, q_pos, body.upper(), ''))
            r_pos += len(body)
        if kind != '*' and len(body) <= _LARGE_INDEL:
            small_diff_bases += len(body)

    def _show(bases: str) -> str:
        return bases if len(bases) <= 20 else f'{len(bases):,} bp'

    aligned_ref = int(hit[8]) - int(hit[7])
    return {
        'ref_len': len(ref_backbone),
        'ref_covered': aligned_ref / len(ref_backbone),
        # identity over the aligned bases, leaving out no-calls and large indels
        'identity': matches / (matches + small_diff_bases),
        'n_substitutions': sum(v[0] == 'Substitution' for v in variants),
        'n_small_indels': sum(v[0] != 'Substitution' and len(v[3] + v[4]) <= _LARGE_INDEL
                              for v in variants),
        'n_large_indels': sum(len(v[3] + v[4]) > _LARGE_INDEL for v in variants),
        'large_indel': _LARGE_INDEL,
        'variants': [
            {
                'type': kind,
                'size': 1 if kind == 'Substitution' else len(ref_bases + alt),
                # 1-based positions; plasmid position is on the input reference
                'backbone_pos': r + 1,
                'plasmid_pos': (plasmid_offset + r) % plasmid_len + 1,
                'consensus_pos': q + 1,
                'ref': _show(ref_bases) or '-',
                'alt': _show(alt) or '-',
            }
            for kind, r, q, ref_bases, alt in variants
        ],
    }


def backbone_consensus(
    outpath: str,
    backbone_reads_path: str,
    flanks_path: str,
    insert_len: int,
    ref_seq: str,
) -> Optional[dict]:
    """Build a consensus of the plasmid backbone (everything outside the MCS, flanks included).

    Backbones are cut from each anchored target plasmid read by barcode_aligner. Only
    full-length backbones (near the most common backbone length) are used. The starting
    draft is a backbone at that length whose read junction (adapter) is furthest from the
    backbone ends, where polishing could not remove it. The draft is polished by repeatedly
    aligning all backbones to it with minimap2 and calling a new consensus with samtools
    consensus (using the read base qualities), until it stops changing. The consensus starts
    at the downstream (right) flank and runs around the plasmid to the end of the upstream
    (left) flank.

    Writes {exp}.backbone_consensus.fa and {exp}.backbone.bam (backbones aligned to the
    consensus). Returns a dict of results for the report, or None if it could not be built.
    """
    exp_name = os.path.basename(outpath)
    log_path = f'{outpath}/Log.txt'

    reads = _read_fasta(backbone_reads_path)
    os.remove(backbone_reads_path)
    if len(reads) < _MIN_READS:
        print(f'Warning: only {len(reads)} anchored target reads; skipping backbone consensus.')
        return None

    peak_len = _peak_length(np.array([len(seq) for _, seq in reads]))
    kept = [(name, seq) for name, seq in reads
            if abs(len(seq) - peak_len) <= _LENGTH_TOLERANCE * peak_len]
    if len(kept) < _MIN_READS:
        print(f'Warning: only {len(kept)} full-length backbone reads; '
              'skipping backbone consensus.')
        return None

    kept = _with_qualities(kept, f'{outpath}/{exp_name}.fastq')
    if len(kept) < _MIN_READS:
        print('Warning: could not find base qualities for the backbone reads; '
              'skipping backbone consensus.')
        return None

    kept_path = f'{outpath}/.tmp.backbone_reads.fq'
    draft_path = f'{outpath}/.tmp.backbone_draft.fa'
    tmp_bam = f'{outpath}/.tmp.backbone.bam'
    with open(kept_path, 'w') as fh:
        for name, seq, qual, _ in kept:
            fh.write(f'@{name}\n{seq}\n+\n{qual}\n')
    near_peak = [r for r in kept if abs(len(r[1]) - peak_len) <= _PEAK_WINDOW * peak_len]
    name, seq, _, _ = max(near_peak or kept, key=lambda r: (r[3], -abs(len(r[1]) - peak_len)))
    _write_fasta(draft_path, [(name, seq)])

    consensus_path = f'{outpath}/{exp_name}.backbone_consensus.fa'
    for _ in range(_MAX_POLISH_ROUNDS):
        _align_sorted(draft_path, kept_path, tmp_bam, log_path)
        # -a keeps the full draft length (uncovered bases become N); deletions are dropped
        # and insertions kept so the consensus follows the reads rather than the draft
        _run(f'samtools consensus -a --show-del no --show-ins yes -f fasta '
             f'-o {draft_path} {tmp_bam}', log_path)
        previous, seq = seq, _read_fasta(draft_path)[0][1]
        if seq == previous:
            break
    os.remove(tmp_bam)
    os.remove(f'{tmp_bam}.bai')
    os.remove(draft_path)
    header = f'{exp_name}_backbone_consensus length={len(seq)} reads={len(kept)}'
    _write_fasta(consensus_path, [(header, seq)])

    # final alignment of the backbones to the consensus, kept for inspection (e.g. IGV)
    final_bam = f'{outpath}/{exp_name}.backbone.bam'
    _align_sorted(consensus_path, kept_path, final_bam, log_path)
    os.remove(kept_path)
    depth = _depth(final_bam, len(seq))

    ref_backbone = _reference_backbone(ref_seq, flanks_path, insert_len)
    comparison = None
    if ref_backbone is not None:
        comparison = _compare_to_reference(consensus_path, ref_backbone[0], ref_backbone[1],
                                           len(ref_seq), outpath, log_path)

    left, right, _ = _load_flanks(flanks_path, insert_len)
    print(f'Backbone consensus: {len(seq)} bp from {len(kept)} reads\n')
    return {
        'sequence': seq,
        'header': header,
        'fasta_name': os.path.basename(consensus_path),
        'reads_anchored': len(reads),
        'reads_used': len(kept),
        'peak_len': peak_len,
        'length_tolerance': int(100 * _LENGTH_TOLERANCE),
        'mean_depth': float(depth.mean()),
        'min_depth': int(depth.min()),
        'n_count': seq.count('N'),
        'comparison': comparison,
        'depth': depth,
        'left_flank_len': len(left),
        'right_flank_len': len(right),
    }
