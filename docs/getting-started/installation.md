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

!!! warning "Apple Silicon / ARM Linux"
    `parasail` publishes no prebuilt wheel for arm64, so pip will try to compile it from
    source (which needs autoconf, automake, libtool and m4). Installing it from bioconda
    first avoids the build entirely:

    ```bash
    conda install -c bioconda parasail-python
    pip install LongBarcodeQC
    ```

To install the development version directly from GitHub:

```bash
pip install git+https://github.com/goodez/LongBarcodeQC.git
```

After installation, the `lbqc` command will be available in your environment.
