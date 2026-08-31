# Installation

## Recommended: conda

The simplest way to install LongBarcodeQC is with conda. This pulls in the external tools
(minimap2, samtools) and the compiled dependencies (parasail, cutadapt) automatically, on
every supported platform — including Apple Silicon (osx-arm64) and ARM Linux, where
`parasail` has no wheel on PyPI:

```bash
conda install -c conda-forge -c bioconda longbarcodeqc
```

Both channels are required: `longbarcodeqc` and the bioinformatics tools come from
bioconda, while Python and the general-purpose libraries come from conda-forge.

## Alternative: PyPI

LongBarcodeQC is also published on PyPI. Note that `pip` will **not** install minimap2 or
samtools, and `parasail` has no prebuilt arm64 wheel, so install those dependencies with
conda first:

```bash
conda install -c conda-forge -c bioconda minimap2 samtools parasail-python cutadapt
pip install LongBarcodeQC
```

## Development version

```bash
pip install git+https://github.com/ArpiarSaundersLab/LongBarcodeQC.git
```

## Verifying the installation

```bash
lbqc --version
lbqc --help
```
