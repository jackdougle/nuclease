mod kmer;
mod proto_duk;
// mod test;

use clap::Parser;

/// K-mer matching tool for sequence analysis
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// K-mer size
    #[arg(short, long, default_value_t = 21)]
    k: usize,

    /// Threshold for k-mer matches
    #[arg(short, long, default_value_t = 3)]
    threshold: usize,

    /// Use canonical k-mers
    #[arg(short, long, default_value_t = true)]
    canonical: bool,

    /// Reference file path
    #[arg(short, long, default_value_t = String::from("in/references.fa"))]
    reference: String,

    /// Reads file path
    #[arg(short, long, default_value_t = String::from("in/100k_150.fq"))]
    query: String,

    #[arg(short, long, default_value_t = String::from("in/reference_kmers.bin"))]
    serialized_kmers_filename: String,

    #[arg(short, long, default_value_t = String::from("out/matched.fa"))]
    matched_path: String,

    #[arg(short, long, default_value_t = String::from("out/unmatched.fa"))]
    unmatched_path: String,
}

fn main() {
    let args = Args::parse();
    proto_duk::run(args);
}
