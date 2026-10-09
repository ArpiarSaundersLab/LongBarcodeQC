import bisect
import gzip
import json
import os
import random
import re
import subprocess
from importlib.resources import files
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
# at most this many full-length backbones are used, as a reproducible random subsample.
# samtools consensus time and memory grow faster than linearly with depth, while the
# consensus stops changing well before this (500-1,000 reads matched all ~2,000 in tests)
_MAX_READS = 5000
_SAMPLE_SEED = 0
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
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
# read support for a difference compares each read to the reference with and without
# the difference, this far either side of it; the read's bases are taken with a little extra
# on each side (free end gaps) so an indel near the window edge does not decide the vote
_SUPPORT_WINDOW = 10
_SUPPORT_PAD = 6
_SUPPORT_MATRIX = parasail.matrix_create('ACGTN', 2, -3)
# sequence context is read from the reference this far either side of a difference
_CONTEXT_MARGIN = 5
_HOMOPOLYMER_MIN = 5


def _iter_fasta(path: str):
    """Yield (name, sequence) from a FASTA file one record at a time."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, ''.join(chunks).upper()


def _read_fasta(path: str) -> list[tuple[str, str]]:
    return list(_iter_fasta(path))


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
    with gzip.open(fastq_path, 'rt') as fh:
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


def load_reference(name: str) -> dict:
    """Load a bundled backbone reference (AP-Amp, AP-Kan or EV).

    These are cut from the annotated plasmid maps: everything between the A and B transfer
    sites (not including them), oriented like the consensus (right flank first), with the
    map's high-level features in backbone coordinates.
    """
    path = files('longbarcodeqc.plasmids').joinpath(f'backbones/{name}.json')
    return json.loads(path.read_text(encoding='utf-8'))


def reference_from_plasmid(ref_seq: str, flanks_path: str, insert_len: int,
                           name: str) -> Optional[dict]:
    """Build a backbone reference (without features) from a user-provided plasmid.

    Returns None if the plasmid does not contain both flanks.
    """
    # the -p help suggests giving the plasmid concatenated to itself; use one copy
    half = len(ref_seq) // 2
    if len(ref_seq) % 2 == 0 and ref_seq[:half] == ref_seq[half:]:
        ref_seq = ref_seq[:half]
    found = _reference_backbone(ref_seq, flanks_path, insert_len)
    if found is None:
        return None
    backbone, offset = found
    return {'name': name, 'label': os.path.splitext(name)[0], 'source': name,
            'map_length': len(ref_seq), 'map_start': offset + 1,
            'sequence': backbone, 'features': []}


def _reference_backbone(ref_seq: str, flanks_path: str, insert_len: int) -> Optional[tuple[str, int]]:
    """Cut the backbone out of the reference plasmid with the same flank anchoring used on reads.

    Returns (backbone, 0-based plasmid position of the backbone's first base), or None if
    the reference does not contain both flanks.
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


def _context(ref_seq: str, start: int, end: int) -> list[str]:
    """Label sequence that nanopore basecalling often miscalls near a difference.

    Looks at reference bases start..end (0-based, end exclusive; start == end for an
    insertion) plus _CONTEXT_MARGIN either side for Dam (GATC) and Dcm (CCWGG) methylation
    sites, which plasmids grown in dam+/dcm+ E. coli carry, and for homopolymers of at
    least _HOMOPOLYMER_MIN bases (labelled with their full length, e.g. G×6). A site or
    homopolymer counts if any of it falls in that window.
    """
    lo, hi = max(0, start - _CONTEXT_MARGIN), min(len(ref_seq), end + _CONTEXT_MARGIN)
    labels = []
    for label, motif in (('Dam (GATC)', 'GATC'), ('Dcm (CCWGG)', 'CC[AT]GG')):
        # widen the search so a site partly inside the window is found
        offset = max(0, lo - 4)
        if any(offset + m.start() < hi and offset + m.end() > lo
               for m in re.finditer(motif, ref_seq[offset:hi + 4])):
            labels.append(label)
    i = lo
    while i < hi:
        run_start, run_end = i, i
        while run_start > 0 and ref_seq[run_start - 1] == ref_seq[i]:
            run_start -= 1
        while run_end < len(ref_seq) and ref_seq[run_end] == ref_seq[i]:
            run_end += 1
        if run_end - run_start >= _HOMOPOLYMER_MIN:
            labels.append(f'{ref_seq[i]}×{run_end - run_start}')
        i = run_end
    return labels


def _read_support(bam_path: str, cons_len: int, ref_seq: str,
                  diffs: list[tuple[str, int, int, str, str]]) -> list[Optional[float]]:
    """Fraction of the reads covering each difference that carry it.

    For each difference (kind, reference pos, consensus pos, reference bases, consensus
    bases; 0-based), the reference from _SUPPORT_WINDOW before it to _SUPPORT_WINDOW after
    it is compared with and without the difference, so other differences and uncalled
    bases nearby do not sway the vote. Both versions are aligned to the read's bases over
    that stretch (reads are aligned to the consensus), and the read counts for the
    consensus only if it scores higher. Ties stay in the denominator. Comparing whole
    stretches rather than single columns keeps reads that are each wrong in their own way
    (as at methylation sites) from adding up to a call none of them supports. Returns None
    for a difference no read covers.
    """
    sites = []
    for _, r, q, ref_bases, alt in diffs:
        before = ref_seq[max(0, r - _SUPPORT_WINDOW):r]
        after = ref_seq[r + len(ref_bases):r + len(ref_bases) + _SUPPORT_WINDOW]
        lo, hi = max(0, q - len(before)), min(cons_len, q + len(alt) + len(after))
        sites.append((lo, hi, before + alt + after, before + ref_bases + after))
    covered = [0] * len(sites)
    for_consensus = [0] * len(sites)

    proc = subprocess.Popen(['samtools', 'view', '-F', '0x904', bam_path],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        fields = line.split('\t', 11)
        start, cigar, seq = int(fields[3]) - 1, fields[5], fields[9]
        # read position at each consensus position the read spans (a base deleted in the
        # read maps to where the deletion sits), plus one past the aligned end
        ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
        span = sum(n for n, op in ops if op in 'MDN=X')
        read_pos = np.empty(span + 1, dtype=np.int64)
        c, q = 0, 0
        for n, op in ops:
            if op in 'M=X':
                read_pos[c:c + n] = np.arange(q, q + n)
                c += n
                q += n
            elif op in 'DN':
                read_pos[c:c + n] = q
                c += n
            elif op in 'IS':
                if c == span:  # trailing soft clip: keep the end at the last aligned base
                    break
                q += n
        read_pos[span] = q
        end = start + span
        for i, (lo, hi, cons_hap, ref_hap) in enumerate(sites):
            if start > lo or end < hi:
                continue
            segment = seq[read_pos[max(start, lo - _SUPPORT_PAD) - start]:
                          read_pos[min(end, hi + _SUPPORT_PAD) - start]]
            covered[i] += 1
            if not segment:
                continue
            cons_score = parasail.sg_dx(cons_hap, segment, 4, 1, _SUPPORT_MATRIX).score
            ref_score = parasail.sg_dx(ref_hap, segment, 4, 1, _SUPPORT_MATRIX).score
            for_consensus[i] += cons_score > ref_score
    if proc.wait() != 0:
        raise subprocess.CalledProcessError(proc.returncode, f'samtools view {bam_path}')
    return [k / n if n else None for k, n in zip(for_consensus, covered)]


def _compare_to_reference(consensus_path: str, bam_path: str, reference: dict, outpath: str,
                          log_path: str) -> Optional[dict]:
    """Align the consensus to the reference backbone, list the differences with their read
    support and sequence context, and carry the reference features over to consensus
    coordinates."""
    ref_backbone = reference['sequence']
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
    blocks = []  # aligned (reference start, consensus start, length), for lifting positions
    matches = 0
    small_diff_bases = 0
    for op in _CS_OP_RE.findall(cs):
        kind, body = op[0], op[1:]
        if kind == ':':
            blocks.append((r_pos, q_pos, int(body)))
            matches += int(body)
            q_pos += int(body)
            r_pos += int(body)
            continue
        if kind == '*':
            blocks.append((r_pos, q_pos, 1))
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

    q_len = int(hit[1])
    block_starts = [b[0] for b in blocks]

    def lift(r: int) -> int:
        """0-based reference position -> 0-based consensus position."""
        i = bisect.bisect_right(block_starts, r) - 1
        if i < 0:  # before the alignment: extend back from its start
            return max(0, blocks[0][1] - (blocks[0][0] - r))
        b_r, b_q, n = blocks[i]
        if r < b_r + n:
            return b_q + (r - b_r)
        if i + 1 < len(blocks):  # deleted in the consensus: next aligned base
            return blocks[i + 1][1]
        return min(q_len - 1, b_q + n - 1 + (r - (b_r + n - 1)))

    features = [dict(f, start=lift(f['start'] - 1) + 1, end=lift(f['end'] - 1) + 1)
                for f in reference['features']]

    def in_features(r: int) -> str:
        return ', '.join(f['name'] for f in reference['features'] if f['start'] <= r + 1 <= f['end'])

    support = _read_support(bam_path, q_len, ref_backbone, variants)
    aligned_ref = int(hit[8]) - int(hit[7])
    return {
        'ref_label': reference['label'],
        'ref_source': reference['source'],
        'features': features,
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
                # 1-based; the reference backbone and the consensus both start at the right flank
                'backbone_pos': r + 1,
                'feature': in_features(r),
                'consensus_pos': q + 1,
                'ref': _show(ref_bases) or '-',
                'alt': _show(alt) or '-',
                'support': s,
                # large differences span many motifs, so context is only given for small ones
                'context': (_context(ref_backbone, r, r + len(ref_bases))
                            if len(ref_bases + alt) <= _LARGE_INDEL else []),
            }
            for (kind, r, q, ref_bases, alt), s in zip(variants, support)
        ],
    }


def backbone_consensus(
    outpath: str,
    backbone_reads_path: str,
    flanks_path: str,
    insert_len: int,
    reference: Optional[dict],
) -> Optional[dict]:
    """Build a consensus of the plasmid backbone (everything outside the MCS, flanks included).

    Backbones are cut from each anchored target plasmid read by barcode_aligner. Only
    full-length backbones (near the most common backbone length) are used, at most
    _MAX_READS of them (a fixed-seed random subsample, so reruns give the same result). The
    backbone file is read twice (lengths, then the chosen sequences) rather than held in
    memory, since large libraries can have many thousands of 10+ kb backbones. The starting
    draft is a backbone at that length whose read junction (adapter) is furthest from the
    backbone ends, where polishing could not remove it. The draft is polished by repeatedly
    aligning all backbones to it with minimap2 and calling a new consensus with samtools
    consensus (using the read base qualities), until it stops changing. The consensus starts
    at the downstream (right) flank and runs around the plasmid to the end of the upstream
    (left) flank.

    If a reference backbone is given (see load_reference / reference_from_plasmid), the
    consensus is compared to it and its features are carried over to the consensus. Each
    difference is flagged with its read support and sequence context, not corrected.

    Writes {exp}.backbone_consensus.fa and {exp}.backbone.bam (backbones aligned to the
    consensus). Returns a dict of results for the report, or None if it could not be built.
    """
    exp_name = os.path.basename(outpath)
    log_path = f'{outpath}/Log.txt'

    lengths = [(name, len(seq)) for name, seq in _iter_fasta(backbone_reads_path)]
    if len(lengths) < _MIN_READS:
        os.remove(backbone_reads_path)
        print(f'Warning: only {len(lengths)} reads with the MCS found; skipping backbone consensus.')
        return None

    peak_len = _peak_length(np.array([n for _, n in lengths]))
    full_length = [name for name, n in lengths if abs(n - peak_len) <= _LENGTH_TOLERANCE * peak_len]
    if len(full_length) < _MIN_READS:
        os.remove(backbone_reads_path)
        print(f'Warning: only {len(full_length)} full-length backbone reads; '
              'skipping backbone consensus.')
        return None
    chosen = set(full_length)
    if len(full_length) > _MAX_READS:
        # sample from input order (names are <run>_<read number>), not the order reads came out
        # of the sorted alignment, which can vary between runs
        in_order = sorted(full_length, key=lambda name: (len(name), name))
        chosen = set(random.Random(_SAMPLE_SEED).sample(in_order, _MAX_READS))
        print(f'Using a random {_MAX_READS:,} of {len(full_length):,} full-length backbones')
    kept = [(name, seq) for name, seq in _iter_fasta(backbone_reads_path) if name in chosen]
    os.remove(backbone_reads_path)

    kept = _with_qualities(kept, f'{outpath}/{exp_name}.fastq.gz')
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

    comparison = None
    if reference is not None:
        comparison = _compare_to_reference(consensus_path, final_bam, reference, outpath, log_path)

    left, right, _ = _load_flanks(flanks_path, insert_len)
    print(f'Backbone consensus: {len(seq)} bp from {len(kept)} reads\n')
    return {
        'sequence': seq,
        'header': header,
        'fasta_name': os.path.basename(consensus_path),
        'reads_anchored': len(lengths),
        'reads_full_length': len(full_length),
        'reads_used': len(kept),
        'max_reads': _MAX_READS,
        'peak_len': peak_len,
        'length_tolerance': int(100 * _LENGTH_TOLERANCE),
        'mean_depth': float(depth.mean()),
        'min_depth': int(depth.min()),
        'n_count': seq.count('N'),
        'comparison': comparison,
        'features': comparison['features'] if comparison else [],
        'depth': depth,
        'left_flank_len': len(left),
        'right_flank_len': len(right),
    }
