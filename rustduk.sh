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

#!/usr/bin/env bash
set -euo pipefail

# Default max memory usage (85% of system memory)
DEFAULT_MEM_PERCENT=85

# Find the value following '--maxmem' in the argument list
MAX_MEM=""
for ((i=1; i<=$#; i++)); do
    if [[ "${!i}" == "--maxmem" ]]; then
        next=$((i+1))
        MAX_MEM="${!next}"
        break
    fi
done

# Detect available system memory in kilobytes
if [[ "$(uname)" == "Linux" ]]; then
    AVAILABLE_MEM_KB=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
    AVAILABLE_MEM_BYTES=$((AVAILABLE_MEM_KB * 1024))
elif [[ "$(uname)" == "Darwin" ]]; then
    # macOS: get free memory in bytes
    AVAILABLE_MEM_BYTES=$(vm_stat | awk '
        /Pages free/ {free=$3}
        /Pages inactive/ {inactive=$3}
        END {print (free+inactive)*4096}
    ')
else
    echo "Unsupported OS for automatic memory detection."
    AVAILABLE_MEM_BYTES=0
fi

# If user provided memory, parse it. Otherwise, use 85% of available system memory
if [[ -z "$MAX_MEM" ]]; then
    LIMIT_BYTES=$((AVAILABLE_MEM_BYTES * DEFAULT_MEM_PERCENT / 100))
else
    # Convert human-readable input like "4G" or "800M" to bytes
    case "$MAX_MEM" in
        *G) LIMIT_BYTES=$(( ${MAX_MEM%G} * 1024 * 1024 * 1024 )) ;;
        *M) LIMIT_BYTES=$(( ${MAX_MEM%M} * 1024 * 1024 )) ;;
        *K) LIMIT_BYTES=$(( ${MAX_MEM%K} * 1024 )) ;;
        *) LIMIT_BYTES=$MAX_MEM ;;  # Assume raw bytes
    esac
fi

echo "Limiting memory to $LIMIT_BYTES bytes."

# Set the memory limit for this shell and any processes it spawns
if [[ "$(uname)" == "Linux" ]]; then
    ulimit -v $((LIMIT_BYTES / 1024))
fi

echo "Running: rust-duk $@" >&2

./target/release/rust-duk "$@"