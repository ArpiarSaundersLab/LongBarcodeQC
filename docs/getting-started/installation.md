# Installation

LongBarcodeQC requires **Python ≥ 3.10** and several external bioinformatics tools
(minimap2, samtools, parasail, cutadapt). These are installed with conda, which provides
prebuilt binaries for every supported platform — including Apple Silicon (osx-arm64) and
ARM Linux, where `parasail` and `cutadapt` have no wheels on PyPI and would otherwise be
compiled from source.

## 1. Create the environment

```bash
conda create -n longbarcodeqc -c conda-forge -c bioconda \
  python=3.12 minimap2 samtools parasail-python cutadapt
conda activate longbarcodeqc
```

Or use the environment file from the repository:

```bash
conda env create -f environment.yml
conda activate longbarcodeqc
```

## 2. Install LongBarcodeQC

```bash
pip install LongBarcodeQC
```

To install the development version instead of the released one:

```bash
pip install git+https://github.com/ArpiarSaundersLab/LongBarcodeQC.git
```

After installation, the `lbqc` command will be available in your environment.

## Verifying the installation

```bash
lbqc --help
```
