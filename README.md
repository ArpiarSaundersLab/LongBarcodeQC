<img src="https://raw.githubusercontent.com/ArpiarSaundersLab/LongBarcodeQC/main/icon.png" alt="LongBarcodeQC" width="400"/>

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

LongBarcodeQC requires **Python ≥ 3.10** and several external bioinformatics tools
(minimap2, samtools, parasail, cutadapt). Install them with conda, which provides
prebuilt binaries for all supported platforms including Apple Silicon and ARM Linux:

```bash
conda create -n longbarcodeqc -c conda-forge -c bioconda \
  python=3.12 minimap2 samtools parasail-python cutadapt
conda activate longbarcodeqc
pip install LongBarcodeQC
```

Or create the environment from the file in the repository:

```bash
conda env create -f environment.yml
conda activate longbarcodeqc
pip install LongBarcodeQC
```

To install the development version instead of the released one:

```bash
pip install git+https://github.com/ArpiarSaundersLab/LongBarcodeQC.git
```

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
| `EV` | Expression Vector — 3 sites x 256 barcodes (768 total) |
| `AP` | Assembly Plasmid — 3 sites x 256 barcodes (768 total) |
| `TS` | TritSeq — 4 sites x 4 positions x 3 barcodes (48 total) |

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
| `-T`, `--trim` | Trim the ONT Rapid (RAP) adapter and its leader sequence with cutadapt before alignment |
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
| `<name>.cutadapt.txt` | cutadapt adapter trimming report (only with `-T`) |

## Requirements

- Python ≥ 3.10
- minimap2 and samtools (external tools, installed with conda — see above)
- parasail and cutadapt ≥ 5.2 (compiled dependencies; installed with conda — see above)
- pandas, matplotlib, seaborn, jinja2, pyarrow, tqdm (pure-Python; installed automatically with pip)

`cutadapt` is only needed for the `-T`/`--trim` option, but it is installed as a
dependency so trimming works out of the box.

## Test data

Four small test datasets (3,000 reads each) are included in the GitHub repository under
`longbarcodeqc/test/data/` to verify an installation. They are not shipped in the PyPI
package to keep the download small — clone the repository to use them:

```bash
git clone https://github.com/ArpiarSaundersLab/LongBarcodeQC.git
cd LongBarcodeQC

lbqc -i longbarcodeqc/test/data/PadlockSeq_AP/ -o /tmp/test_AP -b AP -l 300
lbqc -i longbarcodeqc/test/data/PadlockSeq_EV/ -o /tmp/test_EV -b EV -S -a -l 300
lbqc -i longbarcodeqc/test/data/TritSeq_AP/   -o /tmp/test_TS_AP -b TS -l 300
lbqc -i longbarcodeqc/test/data/TritSeq_EV/   -o /tmp/test_TS_EV -b TS -S -a -l 300
```

Each run writes a `*_summary_report.html` file that can be opened in a browser.
Note that the output directory must be empty.

## Citation

If you use LongBarcodeQC in your work, please cite:

> Goode Z, et al. LongBarcodeQC. (manuscript in preparation)

## License

MIT — see [LICENSE](LICENSE).
