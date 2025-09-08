mod duk;
mod kmer_processor;

use clap::Parser;
#[cfg(unix)]
use libc::{RLIM_INFINITY, RLIMIT_AS, RLIMIT_RSS, getrlimit, rlimit, setrlimit};
use std::io;

const HELP: &str = "Nuclease 1.0.0
Written by Jack Douglass
Last modified September 7, 2025

Purpose: Compares DNA sequences from input file to DNA sequences from reference file
    using k-mer analysis. Splits up reference file sequences into k-mers of specified
    length to build k-mer index, then compares each input read sequence for matching
    k-mers. If a read sequence has >= minhits matching k-mers, it will be printed as a
    match. By virtue of being in Rust, very memory-efficient and performant. Processes
    paired reads, in two files or as a single interleaved file.

Input Parameters
    --in <file>         Input FASTA/FASTQ file containing reads to be filtered. Use
                        'stdin.fq' or 'stdin' to pipe from stdin.
    --in2 <file>        (Optional) Second input file for 2nd pair of reads. Must be
                        same length as main input file.
    --ref <file>        Reference FASTA/FASTQ file containing sequences to build
                        reference k-mer index. Program will serialize reference
                        k-mers and build to --binref path for future use. Not
                        necessary if '--binref <file>' is provided.
    --binref <file>     (Optional) Binary file containing serialized k-mer index,
                        increases performance. Nuclease makes this automatically based
                        on '--ref' file if path is also given here. Increases speed.

Output Parameters: use 'stdout.fa / stdout.fq to pipe to stdout'
    --outm <file>       Output file for reads that have >= minhits k-mer matches to
                        reference k-mers. FASTA format if .fa, .fna, or .fasta, FASTQ
                        format otherwise.
    --outu <file>       Output file for reads that have < minhits k-mer matches to
                        reference k-mers. FASTA format if .fa, .fna, or .fasta, FASTQ
                        format otherwise.
    --outm2 <file>      (Optional) Output file for 2nd pair of matched reads.
    --outu2 <file>      (Optional) Output file for 2nd pair of unmatched reads.

Memory & Performance Parameters
    --k 21              K-mer size (number of bases per k-mer). Ranges from 1-31,
                        larger k-mers will have less matches.
    --minhits 1         Minimum number of k-mer matches a read sequence must have to
                        be considered a match.
    --threads auto      Number of threads to use for parallel processing.
    --maxmem auto       Maximum memory to use. Program will use ~50% of available
                        memory by default. '--maxmem 5G' will specify 5 gigabytes,
                        '--maxmem 200M' will specify 200 megabytes. Memory limiting is
                        only available on Linux.
    --interleaved       Enable flag to input as interleaved paired-end reads, omit flag
                        for unpaired reads.

Function and usage documentation at ./README.md
Contact jack.gdouglass@gmail.com for any questions or issues encountered.
";

#[derive(Parser)]
#[command(version, before_help = "A fast DNA decontamination Rust program using k-mers", long_about = HELP)]
struct Args {
    /// Amount of bases in a k-mer
    #[arg(long)]
    k: Option<usize>,

    /// Min number of k-mer hits to match a read
    #[arg(long)]
    minhits: Option<u8>,

    /// Max amount of threads to use
    #[arg(long)]
    threads: Option<usize>,

    /// Memory cap in human-readable format
    #[arg(long)]
    maxmem: Option<String>,

    /// FASTA/FASTQ path for reference sequences
    #[arg(long)]
    r#ref: String,

    /// FASTA/FASTQ path for read sequences
    #[arg(long)]
    r#in: String,

    /// FASTA/FASTQ path for 2nd pair of reads
    #[arg(long)]
    in2: Option<String>,

    /// Binary file containing serialized ref k-mers
    #[arg(long)]
    binref: Option<String>,

    /// Output file of matched reads
    #[arg(long)]
    outm: String,

    /// Output file of unmatched reads
    #[arg(long)]
    outu: String,

    /// Output file for second pair of matched reads
    #[arg(long)]
    outm2: Option<String>,

    /// Output file for second pair of unmatched reads
    #[arg(long)]
    outu2: Option<String>,

    /// Enabling flag signals interleaved input
    #[arg(long)]
    interinput: bool,
}

fn parse_memory_size(size_str: &str) -> Result<u64, String> {
    let size_str = size_str.trim().to_uppercase();

    if let Some(num_str) = size_str.strip_suffix('G') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024 * 1024 * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else if let Some(num_str) = size_str.strip_suffix('M') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024 * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else if let Some(num_str) = size_str.strip_suffix('K') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else {
        size_str
            .parse::<u64>()
            .map_err(|_| format!("Invalid size format: {}", size_str))
    }
}

#[cfg(unix)]
fn try_set_memory_limit(max_bytes: u64) -> io::Result<()> {
    // First try RLIMIT_RSS (physical memory) - more reliable on macOS
    if try_set_limit(RLIMIT_RSS, max_bytes, "physical").is_ok() {
        return Ok(());
    }

    // If RSS fails, try RLIMIT_AS (virtual memory) - works better on Linux
    if try_set_limit(RLIMIT_AS, max_bytes, "virtual").is_ok() {
        return Ok(());
    }

    Err(io::Error::new(
        io::ErrorKind::Unsupported,
        "System does not support rlimit memory limiting",
    ))
}

#[cfg(unix)]
fn try_set_limit(resource: i32, max_bytes: u64, mem_type: &str) -> io::Result<()> {
    // Get current limits
    let mut current_limit = rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };

    let res = unsafe { getrlimit(resource, &mut current_limit) };
    if res != 0 {
        return Err(io::Error::last_os_error());
    }

    if current_limit.rlim_cur != RLIM_INFINITY && current_limit.rlim_cur <= max_bytes {
        println!(
            "Memory limit already set to {:.1} MB via {}",
            current_limit.rlim_cur as f64 / 1024.0 / 1024.0,
            mem_type
        );
        return Ok(());
    }

    // Set the new limit
    let new_limit = rlimit {
        rlim_cur: max_bytes,
        rlim_max: current_limit.rlim_max,
    };

    println!(
        "Attempting to set {} memory limit to {:.1} MB",
        mem_type,
        max_bytes as f64 / 1024.0 / 1024.0
    );

    let res = unsafe { setrlimit(resource, &new_limit) };
    if res != 0 {
        let err = io::Error::last_os_error();
        return Err(err);
    }

    println!(
        "Successfully set {} limit to {:.1} MB",
        mem_type,
        max_bytes as f64 / 1024.0 / 1024.0
    );
    Ok(())
}

#[cfg(not(unix))]
fn try_set_memory_limit(_max_bytes: u64) -> io::Result<()> {
    println!("Memory limiting only works on UNIX systems\nRunning with default memory allocation");
    Ok(()) // Don't fail, just warn
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    if let Some(ref mem_str) = args.maxmem {
        match parse_memory_size(mem_str) {
            Ok(max_memory) => match try_set_memory_limit(max_memory) {
                Ok(()) => {}
                Err(e) => {
                    eprintln!("Could not set memory limit: {}", e);
                    eprintln!("Running with default memory allocation");
                }
            },
            Err(e) => {
                eprintln!("Error parsing memory size '{}': {}", mem_str, e);
                eprintln!(
                    "Valid formats: 'XG', 'XM', 'XK', or number with no suffix representing bytes"
                );
                std::process::exit(1);
            }
        }
    }

    duk::run(args)?;

    Ok(())
}
