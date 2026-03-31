# Installation

## 1. Install external dependencies

LongBarcodeQC requires **minimap2** and **samtools** to be installed separately.

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

## 2. Install LongBarcodeQC

```bash
pip install LongBarcodeQC
```

To install directly from GitHub:

```bash
pip install git+ssh://git@github.com/goodez/LongBarcodeQC.git
```

After installation, the `lbqc` command will be available in your environment.
