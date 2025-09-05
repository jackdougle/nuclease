#!/bin/bash

intro() {
echo "
Written by Jack Douglass
Last modified September 5, 2025

Purpose: 

Usage: rustduk.sh --ref <file> --in <file> --outu <file> --outm <file> --k <int> ...

Input Parameters:
    --in <file>         Input .fa/.fq file containing reads to be filtered.
    --in2 <file>        (Optional) Second input file for 2nd pair of reads.
    --ref <file>        Reference .fa/.fq file containing sequences to build 
                        reference k-mer index. Program will serialize reference
                        k-mers and build to --binref path for future use. Not 
                        necessary if '--binref <file>' is provided.
    --binref <file>     (Optional) Binary file containing serialized k-mer index,
                        increases performance.

Output Parameters:
    --outm <file>       Output file for reads that have >= minhits k-mer matches to
                        reference k-mers.
    --outu <file>       Output file for reads that have < minhits k-mer matches to
                        reference k-mers.
    --outm2 <file>      (Optional) Second output file for matched paired-end reads.
    --outu2 <file>      (Optional) Second output file for unmatched paired-end 
                        reads.

Memory & Performance Parameters:
    --k 21              K-mer size (number of bases per k-mer).
    --minhits 1         Minimum number of matching k-mers required for a read to be
                        considered a match.
    --threads auto      Number of threads to use for parallel processing.
    --maxmem auto       Maximum memory to use. Program will use ~87.5% of available
                        memory by default. '--maxmem 5G' will specify 5 gigabytes,
                        '--maxmem 200M' will specify 200 megabytes.
    --interleaved       Treat input as interleaved paired-end reads, omit flag for
                        unpaired reads.

Function and usage documentation at ./README.md
Contact jack.gdouglass@gmail.com for any questions or issues encountered.
"
}

set -euo pipefail

echo "Running: rust-duk $@" >&2

rust-duk "$@"