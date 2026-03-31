# Quickstart

## Basic usage

```bash
lbqc -i <fastq_pass_dir> -o <output_dir> -b <barcodes>
```

## Using a built-in barcode preset

Three barcode libraries are included. Pass the preset keyword to `-b`:

| Keyword | Description |
|---------|-------------|
| `EV` | Expression Vector (256 barcodes) |
| `AP` | Assembly Plasmid (256 barcodes) |
| `TS` | TritSeq (48 barcodes) |

```bash
lbqc -i fastq_pass/ -o results/ -b EV
```

## Using a custom barcode file

Provide a path to a FASTA file where each entry is one barcode sequence and all headers are unique:

```bash
lbqc -i fastq_pass/ -o results/ -b /path/to/barcodes.fa
```

## Full example

```bash
lbqc \
  -i /data/run01/fastq_pass/ \
  -o /results/run01_EV/ \
  -b EV \
  -r enzymes.txt \
  -S
```

See the [CLI Reference](../usage/cli.md) for all available options.
