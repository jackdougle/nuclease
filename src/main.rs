mod engine;
mod kmer_processor;
#[cfg(test)]
mod test;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// K-mer size
    #[arg(short, long, default_value_t = 0)]
    k: usize,

    /// Threshold for k-mer matches
    #[arg(short, long, default_value_t = 0)]
    minhits: u8,

    /// Amount of threads to use
    #[arg(short, long, default_value_t = 0)]
    threads: u8,

    /// Memory cap
    #[arg(short, long, default_value_t = 0)]
    maxmem: u32,

    /// Reference file path
    #[arg(short, long, default_value_t = String::new())]
    r#ref: String,

    #[arg(short, long, default_value_t = String::new())]
    r#in: String,

    #[arg(short, long, default_value_t = String::new())]
    in2: String,

    #[arg(short, long, default_value_t = String::new())]
    binref: String,

    #[arg(short, long, default_value_t = String::new())]
    outm: String,

    #[arg(short, long, default_value_t = String::new())]
    outu: String,

    #[arg(short, long, default_value_t = String::new())]
    outm2: String,

    #[arg(short, long, default_value_t = String::new())]
    outu2: String,

    #[arg(short, long, default_value_t = false)]
    interleaved: bool,
}

fn main() {
    let args = Args::parse();
    engine::run(args);
}
