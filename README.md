<img src="icon.png" alt="LongBarcodeQC" width="400"/>

# LongBarcodeQC

**LongBarcodeQC** is a command-line tool for analyzing combinatorial barcode libraries sequenced with Oxford Nanopore or other long read platforms. It aligns reads, extracts the multi-cloning site (MCS) region, scores each read against a barcode library, and generates an interactive HTML summary report.

## Features

- Merges and preprocesses FASTQ input (one or multiple files)
- Aligns reads to a reference plasmid using minimap2
- Extracts the MCS region from each read using flanking sequence anchors (parasail Smith-Waterman)
- Scores each read against a barcode library and calls the best-matching barcode per position
- Detects restriction enzyme cut sites within the MCS
- Generates a self-contained HTML report with interactive plots (read length distributions, barcode heatmaps, z-score distributions, restriction site summaries)
- Outputs a compressed summary CSV with per-read barcode calls and QC metrics

## Installation

LongBarcodeQC requires **Python ≥ 3.10** and two external bioinformatics tools that must be installed separately.

Install LongBarcodeQC

```bash
pip install LongBarcodeQC
```

Install external dependencies

**minimap2** — long-read aligner:

```bash
# conda (macOS/Linux)
conda install -c bioconda minimap2

# Homebrew (macOS)
brew install minimap2

# apt (Linux)
sudo apt install minimap2
```

**samtools** — SAM/BAM processing:

```bash
# conda (macOS/Linux)
conda install -c bioconda samtools

# Homebrew (macOS)
brew install samtools

# apt (Linux)
sudo apt install samtools
```



To install directly from GitHub:


After installation, the `lbqc` command will be available in your environment.

## Usage

```
lbqc -i <fastq_pass_dir> -o <output_dir> -b <barcodes>
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Path to Nanopore `fastq_pass` directory (or a single FASTQ file) |
| `-o`, `--output` | Output directory (created if it does not exist; must be empty) |
| `-b`, `--barcodes` | Barcode FASTA file — use a preset keyword or a path to a custom FASTA (see below) |

### Barcode presets

Three built-in barcode libraries are included:

| Keyword | Description |
|---------|-------------|
| `EV` | Expression Vector (256 barcodes per site) |
| `AP` | Assembly Plasmid (256 barcodes per site) |
| `TS` | TritSeq (3 barcodes per site) |

```bash
lbqc -i fastq_pass/ -o results/ -b EV
```

To use a custom barcode library, provide a path to a FASTA file where each entry is one barcode sequence and all headers are unique:

```bash
lbqc -i fastq_pass/ -o results/ -b /path/to/barcodes.fa
```

### Optional arguments

| Argument | Description |
|----------|-------------|
| `-p`, `--plasmid` | Path to a custom plasmid FASTA. Default: built-in Assembly Plasmid |
| `-l`, `--insert_length` | Expected insert size in bp. Default: 300 |
| `-f`, `--flanks` | FASTA with upstream and downstream MCS flanking sequences (required for custom plasmids) |
| `-r`, `--enzymes` | Text file listing desired restriction enzyme names and sequences (one per line, comma-separated) |
| `-a`, `--AP` | Flag Assembly Plasmid reads as contamination (useful after transfer to Expression Vector) |
| `-S`, `--SBARRO` | Use SBARRO mode (rabies genome; inserts NNN sequence into MCS for alignment) |
| `-z`, `--zscore` | Manually set z-score threshold for barcode calling (recommended - check html report after initial run) |
| `-N`, `--expected_insertions` | Expected number of insertions per library member (used in read length histogram) |
| `--full-output` | Write full per-barcode alignment score table as a Parquet file |

### Example

```bash
lbqc \
  -i /data/run01/fastq_pass/ \
  -o /results/run01_EV/ \
  -b EV \
  -r enzymes.txt \
  -S \
  -a
```

## Output

| File | Description |
|------|-------------|
| `<name>_summary.csv.gz` | Per-read barcode calls, MCS metrics, and QC flags |
| `<name>_summary_report.html` | Self-contained interactive HTML report |
| `<name>.aligned.fa.gz` | Reads that aligned to the reference plasmid |
| `<name>.unaligned.fa.gz` | Reads that did not align |
| `<name>.parquet` | Full barcode alignment scores per read (only with `--full-output`) |

## Requirements

- Python ≥ 3.10
- minimap2 (external, see installation above)
- samtools (external, see installation above)
- pandas, parasail, matplotlib, seaborn, jinja2, pyarrow, tqdm (installed automatically with pip)
