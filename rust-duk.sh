#!/bin/bash
# rust-duk.sh --ref <file> --in <file> --outu <file> --outm <file> --k <int> ...

intro() {
echo "
Written by Jack Douglass
Last modified September 1, 2025

Purpose:
    Runs the rust-duk k-mer matching tool.

Input Parameters:
    --in <file>         Input FASTA/FASTQ file containing reads to be filtered.
    --in2 <file>        (Optional) Second input file for paired-end reads.
    --ref <file>        Reference FASTA/FASTQ file containing sequences to build k-mer index.
    --binref <file>     (Optional) Path to binary file for serialized k-mer index (for faster startup).

Output Parameters:
    --outm <file>       Output file for reads that match the reference k-mers.
    --outu <file>       Output file for reads that do NOT match the reference k-mers.
    --outm2 <file>      (Optional) Second output file for matched paired-end reads.
    --outu2 <file>      (Optional) Second output file for unmatched paired-end reads.

Memory & Performance Parameters:
    --k 21              K-mer size (number of bases per k-mer).
    --minhits 1         Minimum number of matching k-mers required for a read to be considered a match.
    --threads <int>     Number of threads to use for parallel processing.
    --maxmem <int>      Maximum memory (in MB) to use.
    --interinput        Treat input as interleaved paired-end reads, omit flag for normal reads.

Function and usage documentation at ./README.md
Please contact jack.gdouglass@gmail.com if you encounter any problems.
"
}

echo "$@"

rust-duk "$@"