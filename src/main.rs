mod duk;
mod kmer_processor;
// #[cfg(test)]
// mod test;

use clap::Parser;
use std::time::Instant;

/// K-mer matching tool for sequence analysis
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// K-mer size
    #[arg(short, long, default_value_t = 21)]
    k: usize,

    /// Threshold for k-mer matches
    #[arg(short, long, default_value_t = 1)]
    num_hits: u8,

    /// Amount of threads to use
    #[arg(short, long, default_value_t = 8)]
    threads: u8,

    /// Reference file path
    #[arg(short, long, default_value_t = String::from("in/10k_150.fa"))]
    ref_path: String,

    /// Reads file path
    #[arg(short, long, default_value_t = String::from("in/1m_150.fq"))]
    in_path: String,

    #[arg(short, long, default_value_t = String::from("in/serialized_kmers.bin"))]
    bin_kmers_path: String,

    #[arg(short, long, default_value_t = String::from("out/matched.fa"))]
    matched_path: String,

    #[arg(short, long, default_value_t = String::from("out/unmatched.fa"))]
    unmatched_path: String,

    #[arg(short, long, default_value_t = false)]
    paired_reads: bool,
}

fn main() {
    let start_time = Instant::now();
    let args = Args::parse();
    duk::run(args, start_time);
}
