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

## Trimming adapters

Nanopore reads carry the ONT Rapid (RAP) adapter behind an inconsistently basecalled leader
sequence, adding roughly 90 bp of non-biological sequence to the front of most reads. Use `-T`
to remove the adapter and its leader with cutadapt before alignment:

```bash
lbqc -i fastq_pass/ -o results/ -b EV -T
```

A `<name>.cutadapt.txt` report is written to the output directory, and the trimming stats are
printed to the console.

## Full example

```bash
lbqc \
  -i /data/run01/fastq_pass/ \
  -o /results/run01_EV/ \
  -b EV \
  -r enzymes.txt \
  -T \
  -S
```

See the [CLI Reference](../usage/cli.md) for all available options.
