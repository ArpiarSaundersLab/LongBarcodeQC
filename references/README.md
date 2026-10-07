# Backbone references

Source maps and build script for the backbone references used by the `-c/--consensus` option.

| Reference | Used for | Source map (`maps/`) |
|-----------|----------|----------------------|
| `AP-Amp` | Default plasmids, when AP-Amp has more reads than AP-Kan | `AP-Amp.dna` |
| `AP-Kan` | Default plasmids, when AP-Kan has more reads than AP-Amp | `AP-Kan.dna` |
| `EV` | `-S` runs (any barcode set) | `CVS-N2c(deltaG)-N2c-nl.EGFP-SypGFP-TermBAv2_TspMI_PadlockSeq_A.dna` |

Each reference is everything between the A and B transfer sites, not including them, oriented
like the consensus: from the right MCS flank around the plasmid to the end of the left flank.
The maps are copies of the annotated plasmid maps in the manuscript's supplementary data (SDF1).

To rebuild after changing a map, run from the repository root:

```bash
python references/build_backbone_references.py
```

This writes `longbarcodeqc/plasmids/backbones/<name>.json` with the backbone sequence, its
high-level features and the source map's SHA-256. The script stops if a backbone does not start
and end on the flank sequences in `longbarcodeqc/plasmids/`.
