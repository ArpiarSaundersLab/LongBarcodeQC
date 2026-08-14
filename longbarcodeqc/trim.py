import os
import re
import subprocess

# The ONT Rapid adapter. It sits ~40-50 bp into the read, behind a leader
# sequence that basecalling renders inconsistently. cutadapt's -g removes
# the adapter and everything before it, so the leader goes too.
RAP = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"


def reverse_complement(seq):
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(pairs[base] for base in reversed(seq.upper()))


def find_number(pattern, text):
    """Pull one number out of the report, e.g. '119,372' -> 119372."""
    match = re.search(pattern, text)
    if match:
        return int(match.group(1).replace(",", ""))
    return 0


def parse_report(report):
    """Print the basic numbers from a cutadapt report and return the counts.

    mean_removed is the average bases cut per input read, which is a rough
    check on the size of the adapter + leader block (~92 bp for RAP).
    """
    reads_in = find_number(r"Total reads processed:\s*([\d,]+)", report)
    with_adapter = find_number(r"Reads with adapters:\s*([\d,]+)", report)
    reads_out = find_number(r"Reads written \(passing filters\):\s*([\d,]+)", report)
    bp_in = find_number(r"Total basepairs processed:\s*([\d,]+)", report)
    bp_out = find_number(r"Total written \(filtered\):\s*([\d,]+)", report)

    if reads_in > 0:
        percent = 100.0 * with_adapter / reads_in
        mean_removed = (bp_in - bp_out) / reads_in
    else:
        percent = 0.0
        mean_removed = 0.0

    print("  reads in:            %d" % reads_in)
    print("  reads with adapter:  %d (%.1f%%)" % (with_adapter, percent))
    print("  reads out:           %d" % reads_out)
    print("  mean bases removed:  %.1f" % mean_removed)

    return {"reads_in": reads_in, "reads_out": reads_out}


def trim_reads(fastq_in, fastq_out, log_file=None, threads=1,
               error_rate=0.25, overlap=20, min_length=1, adapter=RAP):
    """Run cutadapt on one FASTQ file, print the stats, and return the counts.

    Pass fastq_out="/dev/null" to get the stats without writing any reads.
    """
    command = [
        "cutadapt",
        "-g", adapter,
        "-a", reverse_complement(adapter),
        "-e", str(error_rate),
        "-O", str(overlap),
        "--times", "2",
        "-m", str(min_length),
        "-j", str(threads),
        "-o", str(fastq_out),
        str(fastq_in),
    ]

    result = subprocess.run(command, capture_output=True, text=True)
    report = result.stdout + result.stderr

    if result.returncode != 0:
        raise RuntimeError("cutadapt failed on %s:\n%s" % (fastq_in, report))

    if log_file:
        with open(log_file, "w") as handle:
            handle.write(report)

    print("trimmed %s" % fastq_in)
    return parse_report(report)


def trim_fastq_in_place(fastq_path, log_file=None, threads=1):
    """Trim adapters and replace the input FASTQ with the trimmed reads.

    Returns the number of reads dropped during trimming (reads that were
    left too short to keep), so they can be reported alongside the
    alignment counts.
    """
    trimmed_path = fastq_path + ".trimmed"
    counts = trim_reads(fastq_path, trimmed_path, log_file=log_file, threads=threads)
    os.replace(trimmed_path, fastq_path)
    return counts["reads_in"] - counts["reads_out"]
