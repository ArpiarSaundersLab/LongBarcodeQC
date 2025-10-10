# CLI Reference

```text
usage: LongBarcodeQC [-h] -i INPUT -o OUTPUT [-p PLASMID] -b BARCODES [-l INSERT_LENGTH] [-f FLANKS] [-r ENZYMES] [-a] [-S] [-z ZSCORE]
                 [-N EXPECTED_INSERTIONS]

Analyze long-read barcodes and generate a summary report

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input fastq files (nanopore fastq_pass folder).
  -o OUTPUT, --output OUTPUT
                        Path to output directory. An output folder will be created if it does not exist. Program will exit if files already exist here.
  -p PLASMID, --plasmid PLASMID
                        Path to plasmid fasta file. The default plasmid is a.31. If using another plasmid, use this option with the path to a fasta file
                        containing your plasmid. For best results, concatenate the circular plasmid sequence to itself to optimize alignment.
  -b BARCODES, --barcodes BARCODES
                        Path to barcode fasta file. One barcode per line. Fasta headers must be unique names.
  -l INSERT_LENGTH, --insert_length INSERT_LENGTH
                        Integer length of expected insert size cloned into MCS. A rough estimate is sufficient. Default: 300 bp
  -f FLANKS, --flanks FLANKS
                        Path to fasta file containing MCS flanking regions. The first sequence must be the upstream flank, and the second sequence must be
                        the downstream flank. Must be provided in the same file.
  -r ENZYMES, --enzymes ENZYMES
                        Path to text file with restriction site names and their sequences. Define one site per line, with the name and sequence separated
                        by comma.
  -a, --a31             Option to add a31 alignment. Mainly used to detect contamination from a31 plasmid. This will flag a31 reads in the final output
                        and score barcodes against them.
  -S, --SBARRO          Use this option if the plasmid is part of the SBARRO system (derived from c.18). Reference will be generated with NNN sequence
                        inserted into the MCS (length of NNNs equal to provided insert length).
  -z ZSCORE, --zscore ZSCORE
                        Optionally set z-score threshold for barcode calling. The default threshold is set by using the z-score for the worst performing
                        barcode alignment in the library (absolute value). Sometimes this value is either too conservative or too permissive. See the
                        z-score distribution in the html output summary to evaluate this. If desired - use this option to manually set the z-score
                        threshold to be a fixed value (e.g. 2.5).
  -N EXPECTED_INSERTIONS, --expected_insertions EXPECTED_INSERTIONS
                        Optionally set the expected number of insertions in a library. This is used in the html report to categorize reads in the read
                        length histogram. By default, this is set to the maximum number of sites/positions found in the dataset.
```

