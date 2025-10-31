mod core;
mod kmer_ops;

use clap::Parser;
use std::io;
use std::time::Instant;

const ABOUT: &str = "Nucleaze 1.2.0
Written by Jack Douglass
Last modified October 25, 2025

Nucleaze compares DNA sequences from input file to DNA sequences from reference
 file using k-mer analysis. Splits up reference file sequences into k-mers of
 specified length to build k-mer index, then compares each input read sequence
 for matching k-mers. If a read sequence has >= minhits matching k-mers, it
 will be printed as a match. Very memory-efficient and performant. Processes
 paired reads in two files or as a single interleaved file.

USAGE: nucleaze --in <reads file> --ref <ref file> ...

INPUT PARAMETERS
    --in <file>         Input FASTA/FASTQ file containing reads to be filtered.
                        Use 'stdin.fq' or 'stdin' to pipe from stdin.
    --in2 <file>        (Optional) Second input file for 2nd pair of reads. 
                        Must be same length as main input file.
    --ref <file>        (-r) Reference FASTA/FASTQ file containing sequences to
                        build reference k-mer index. Program will serialize
                        reference k-mers and build to --binref path for future
                        use. Not necessary if '--binref <file>' is provided.
    --saveref <file>    (-s) Path at which to store serialized k-mer index if
                        no valid serialized k-mer index is provided.
    --binref <file>     (-b) Binary file containing serialized k-mer index,
                        increases performance. Nucleaze makes this
                        automatically based on '--ref' file if '--saveref'
                        <file> is provided.

OUTPUT PARAMETERS: use 'stdout.fa / stdout.fq to pipe to stdout'
    --outm <file>       Output file for reads that have >= minhits k-mer 
                        matches to reference k-mers. FASTA format if .fa, .fna,
                        or .fasta, FASTQ format otherwise.
    --outu <file>       Output file for reads that have < minhits k-mer matches
                        to reference k-mers. FASTA format if .fa, .fna, or
                        .fasta, FASTQ format otherwise.
    --outm2 <file>      Output file for 2nd pair of matched reads.
    --outu2 <file>      Output file for 2nd pair of unmatched reads.

MEMORY & PERFORMANCE PARAMETERS
    --k 21              (-k) K-mer size (number of bases per k-mer). Ranges
                        from 1-31, larger k-mers will have less matches.
    --minhits 1         Minimum number of k-mer matches a read sequence must
                        have to be considered a match.
    --threads auto      (-t) Number of threads to use for parallel processing.
                        Program will use all available threads by default.
    --maxmem auto       (-m) Maximum memory to use. Program will use ~50% of
                        available memory by default. '--maxmem 5G' will specify
                        5 gigabytes, '--maxmem 200M' will specify 200
                        megabytes. Memory limiting is currently a WIP.
    --interinput        (-i) Enable flag to input as interleaved paired-end
                        reads, omit flag for unpaired reads.
    --order             (-o) Enable flag to get read outputs ordered by
                        sequence ID.

Function and usage documentation at ./README.md
Contact jack.gdouglass@gmail.com for any questions or issues encountered.
";

#[derive(Parser)]
#[command(version, override_help = ABOUT)]
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

    /// Path to store serialized reference index
    #[arg(short, long)]
    saveref: Option<String>,

    /// Binary file containing serialized ref k-mers
    #[arg(long)]
    binref: Option<String>,

    /// Output file of matched reads
    #[arg(long)]
    outm: Option<String>,

    /// Output file of unmatched reads
    #[arg(long)]
    outu: Option<String>,

    /// Output file for second pair of matched reads
    #[arg(long)]
    outm2: Option<String>,

    /// Output file for second pair of unmatched reads
    #[arg(long)]
    outu2: Option<String>,

    /// Enabling flag signals interleaved input
    #[arg(long)]
    interinput: bool,

    /// Enabling flag causes ordered output
    #[arg(short, long)]
    order: bool,
}

fn main() -> io::Result<()> {
    let start_time = Instant::now();

    let args = Args::parse();

    let version = env!("CARGO_PKG_VERSION");
    let user_args: Vec<String> = std::env::args().skip(1).collect();
    println!("Nucleaze {} [{}]", version, user_args.join(" "));

    validate_args(&args)?;

    core::run(args, start_time)?;

    Ok(())
}

fn validate_args(args: &Args) -> io::Result<()> {
    // Check input files aren't the same
    if let Some(ref in2) = args.in2 {
        if args.r#in == *in2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "both inputs (--in and --in2) cannot be the same file",
            ));
        }
    }

    // Check output files aren't the same (if provided)
    if let (Some(outm), Some(outu)) = (&args.outm, &args.outu) {
        if outm == outu {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "matches (--outm) and non-matches (--outu) cannot have the same output path",
            ));
        }
    }

    // Check paired outputs if provided and that they aren't the same
    if let (Some(outm2), Some(outu2)) = (&args.outm2, &args.outu2) {
        if outm2 == outu2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Matches (--outm2) and non-matches (--outu2) cannot have the same output path",
            ));
        }
    }

    Ok(())
}
