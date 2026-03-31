# LongBarcodeQC

**LongBarcodeQC** is a command-line tool for analyzing combinatorial barcode libraries sequenced with Oxford Nanopore long reads. It aligns reads to a reference plasmid, extracts the multi-cloning site (MCS) region, scores each read against a barcode library, and generates an interactive HTML summary report.

## Features

- Merges and preprocesses FASTQ input
- Aligns reads to a reference plasmid using minimap2
- Extracts the MCS region from each read using flanking sequence anchors (parasail Smith-Waterman)
- Scores each read against a barcode library and calls the best-matching barcode per position
- Detects restriction enzyme cut sites within the MCS
- Generates a self-contained HTML report with interactive plots (read length distributions, barcode heatmaps, z-score distributions, restriction site summaries)
- Outputs a compressed summary CSV with per-read barcode calls and QC metrics

## Quick links

- [Installation](getting-started/installation.md)
- [Quickstart](getting-started/quickstart.md)
- [CLI Reference](usage/cli.md)
