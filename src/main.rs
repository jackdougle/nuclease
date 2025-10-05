mod duk;
mod kmer_processor;

use clap::Parser;
use std::io;

const ABOUT: &str = "Nuclease 1.0.3
Written by Jack Douglass
Last modified October 5, 2025

Purpose: Compares DNA sequences from input file to DNA sequences from reference
 file using k-mer analysis. Splits up reference file sequences into k-mers of 
 specified length to build k-mer index, then compares each input read sequence 
 for matching k-mers. If a read sequence has >= minhits matching k-mers, it 
 will be printed as a match. Very memory-efficient and performant. Processes 
 paired reads in two files or as a single interleaved file.

INPUT PARAMETERS
    --in <file>         Input FASTA/FASTQ file containing reads to be filtered.
                        Use 'stdin.fq' or 'stdin' to pipe from stdin.
    --in2 <file>        (Optional) Second input file for 2nd pair of reads. 
                        Must be same length as main input file.
    --ref <file>        Reference FASTA/FASTQ file containing sequences to 
                        build reference k-mer index. Program will serialize
                        reference k-mers and build to --binref path for future
                        use. Not necessary if '--binref <file>' is provided.
    --binref <file>     (Optional) Binary file containing serialized k-mer 
                        index, increases performance. Nuclease makes this
                        automatically based on '--ref' file if path is also
                        given here. Increases speed.

OUTPUT PARAMETERS: use 'stdout.fa / stdout.fq to pipe to stdout'
    --outm <file>       Output file for reads that have >= minhits k-mer 
                        matches to reference k-mers. FASTA format if .fa, .fna,
                        or .fasta, FASTQ format otherwise.
    --outu <file>       Output file for reads that have < minhits k-mer matches
                        to reference k-mers. FASTA format if .fa, .fna, or
                        .fasta, FASTQ format otherwise.
    --outm2 <file>      (Optional) Output file for 2nd pair of matched reads.
    --outu2 <file>      (Optional) Output file for 2nd pair of unmatched reads.

MEMORY & PERFORMANCE PARAMETERS
    --k 21              K-mer size (number of bases per k-mer). Ranges from 
                        1-31, larger k-mers will have less matches.
    --minhits 1         Minimum number of k-mer matches a read sequence must
                        have to be considered a match.
    --threads auto      Number of threads to use for parallel processing.
    --maxmem auto       Maximum memory to use. Program will use ~50% of
                        available memory by default. '--maxmem 5G' will specify
                        5 gigabytes, '--maxmem 200M' will specify 200
                        megabytes. Memory limiting is currently a WIP.
    --interinput        Enable flag to input as interleaved paired-end reads,
                        omit flag for unpaired reads.
    --order             Enable flag to get read outputs ordered by sequence ID.

Function and usage documentation at ./README.md
Contact jack.gdouglass@gmail.com for any questions or issues encountered.
";

#[derive(Parser)]
#[command(
    version,
    before_help = "Fast DNA decontamination Rust program using k-mers",
    override_help = ABOUT,
)]
struct Args {
    /// Amount of bases in a k-mer
    #[arg(short, long)]
    k: Option<usize>,

    /// Max amount of threads to use
    #[arg(short, long)]
    threads: Option<usize>,

    /// Min number of k-mer hits to match a read
    #[arg(long)]
    minhits: Option<u8>,

    /// Memory cap in human-readable format
    #[arg(long)]
    maxmem: Option<String>,

    /// FASTA/FASTQ path for reference sequences
    #[arg(short, long)]
    r#ref: String,

    /// Binary file containing serialized ref k-mers
    #[arg(short, long)]
    binref: Option<String>,

    /// FASTA/FASTQ path for read sequences
    #[arg(long)]
    r#in: String,

    /// FASTA/FASTQ path for 2nd pair of reads
    #[arg(long)]
    in2: Option<String>,

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
    #[arg(short, long)]
    interinput: bool,

    /// Enabling flag causes ordered output
    #[arg(short, long)]
    order: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    duk::run(args)?;

    Ok(())
}
