#!/usr/bin/env python
"""Build the bundled backbone references from the annotated SnapGene maps in maps/.

Each reference is the plasmid backbone: everything between the A and B transfer sites,
not including them. It is written in the same orientation as the backbone consensus,
starting with the right (downstream) MCS flank and running around the plasmid to the end
of the left (upstream) flank, and the script checks that it starts and ends exactly on the
flank sequences LongBarcodeQC uses for that plasmid. High-level map features (genes, CDS,
origins, promoters, enhancers, terminators and other features of at least 100 bp) are kept
in backbone coordinates for the report's plasmid map.

Output: longbarcodeqc/plasmids/backbones/<name>.json, recording the source map and its
SHA-256 so each bundled reference can be traced to an exact file.

Usage (from the repository root): python references/build_backbone_references.py
"""
import hashlib
import json
import os
import struct
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
PLASMIDS = os.path.join(HERE, '..', 'longbarcodeqc', 'plasmids')

# name -> (label shown in the report, source map, flank file the tool uses for this plasmid)
REFERENCES = {
    'AP-Amp': ('AP-Amp', 'AP-Amp.dna', 'AP_flanks.fa'),
    'AP-Kan': ('AP-Kan', 'AP-Kan.dna', 'AP_flanks.fa'),
    'EV': ('CVS-N2c(\u0394G)', 'CVS-N2c(deltaG)-N2c-nl.EGFP-SypGFP-TermBAv2_TspMI_PadlockSeq_A.dna',
           'SBARRO_flanks.fa'),
}

KEEP_TYPES = {'CDS', 'gene', 'rep_origin', 'promoter', 'enhancer', 'terminator', 'misc_feature'}
MIN_FEATURE_LEN = 100
CASSETTE_FEATURES = {'MCS', 'A_Transfer_Cassette', 'B_Transfer_Cassette'}


def read_snapgene(path):
    """Return (sequence, features) from a SnapGene .dna file.

    The file is a series of packets (1-byte type, 4-byte big-endian length, body): type 0
    holds the sequence after one topology byte, type 10 holds the features as XML.
    Feature positions are 1-based and inclusive; a feature that crosses the origin of a
    circular map has start > end.
    """
    data = open(path, 'rb').read()
    i, seq, features = 0, None, []
    while i < len(data):
        ptype = data[i]
        length = struct.unpack('>I', data[i + 1:i + 5])[0]
        body = data[i + 5:i + 5 + length]
        i += 5 + length
        if ptype == 0x00:
            seq = body[1:].decode('ascii').upper()
        elif ptype == 0x0A:
            for f in ET.fromstring(body.decode('utf-8')).iter('Feature'):
                ranges = [[int(x) for x in s.get('range').split('-')]
                          for s in f.iter('Segment') if s.get('range')]
                features.append({
                    'name': f.get('name'),
                    'type': f.get('type'),
                    'strand': {'1': '+', '2': '-'}.get(f.get('directionality'), '.'),
                    'start': ranges[0][0],
                    'end': ranges[-1][1],
                })
    return seq, features


def read_flanks(path):
    records = [line.strip().upper() for line in open(path) if line.strip()]
    return records[1], records[3]  # (left/upstream, right/downstream)


def build(name, label, map_file, flank_file):
    map_path = os.path.join(HERE, 'maps', map_file)
    seq, features = read_snapgene(map_path)
    length = len(seq)

    a = next(f for f in features if f['name'] == 'A_Transfer_Cassette')
    b = next(f for f in features if f['name'] == 'B_Transfer_Cassette')
    first, second = sorted((a, b), key=lambda f: f['start'])  # order differs between AP and EV
    start, end = second['end'] + 1, first['start'] - 1           # map positions of the backbone ends
    bb_len = (end - start) % length + 1
    backbone = (seq * 2)[start - 1:start - 1 + bb_len]

    left, right = read_flanks(os.path.join(PLASMIDS, flank_file))
    assert backbone.startswith(right), f'{name}: backbone does not start with the right flank'
    assert backbone.endswith(left), f'{name}: backbone does not end with the left flank'

    def to_backbone(pos):
        return (pos - start) % length + 1

    kept = []
    for f in features:
        if (f['type'] not in KEEP_TYPES or f['name'] in CASSETTE_FEATURES
                or f['name'].lower().startswith('site')):  # barcode annotations
            continue
        s, e = to_backbone(f['start']), to_backbone(f['end'])
        if s > e or e > bb_len:  # overlaps the cassette
            continue
        kept.append({'name': f['name'], 'type': f['type'], 'strand': f['strand'], 'start': s, 'end': e})

    # join pieces of one feature split at the origin of a map saved as linear (e.g. AmpR)
    kept.sort(key=lambda f: f['start'])
    joined = []
    for f in kept:
        if joined and joined[-1]['name'] == f['name'] and f['start'] <= joined[-1]['end'] + 1:
            joined[-1]['end'] = max(joined[-1]['end'], f['end'])
        else:
            joined.append(dict(f))

    # drop short features and repeat annotations of the same element (keep the longest)
    final = []
    for f in sorted(joined, key=lambda f: f['start'] - f['end']):
        if f['end'] - f['start'] + 1 < MIN_FEATURE_LEN:
            continue
        if any(g['name'] == f['name'] and f['start'] <= g['end'] and g['start'] <= f['end'] for g in final):
            continue
        final.append(f)
    final.sort(key=lambda f: f['start'])

    reference = {
        'name': name,
        'label': label,
        'source': map_file,
        'source_sha256': hashlib.sha256(open(map_path, 'rb').read()).hexdigest(),
        'map_length': length,
        'map_start': start,
        'map_end': end,
        'excluded': f'{first["name"]} to {second["name"]} (map {first["start"]}-{second["end"]})',
        'flanks': flank_file,
        'sequence': backbone,
        'features': final,
    }
    out = os.path.join(PLASMIDS, 'backbones', f'{name}.json')
    with open(out, 'w') as fh:
        json.dump(reference, fh, indent=1, ensure_ascii=False)
        fh.write('\n')
    print(f'{name}: {bb_len:,} bp backbone (map {start}..{end} of {length:,}), '
          f'{len(final)} features -> {os.path.relpath(out)}')


if __name__ == '__main__':
    os.makedirs(os.path.join(PLASMIDS, 'backbones'), exist_ok=True)
    for name, (label, map_file, flank_file) in REFERENCES.items():
        build(name, label, map_file, flank_file)
