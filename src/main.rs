mod duk;
mod kmer_processor;

use clap::Parser;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    /// K-mer size
    #[arg(long, default_value_t = 21)]
    k: usize,

    /// Threshold for k-mer matches
    #[arg(long, default_value_t = 1)]
    minhits: u8,

    /// Amount of threads to use
    #[arg(long)]
    threads: Option<usize>,

    /// Memory cap
    #[arg(long)]
    maxmem: Option<String>,

    /// Reference file path
    #[arg(long, default_value_t = String::new())]
    r#ref: String,

    #[arg(long, default_value_t = String::new())]
    r#in: String,

    #[arg(long, default_value_t = String::new())]
    in2: String,

    #[arg(long, default_value_t = String::new())]
    binref: String,

    #[arg(long, default_value_t = String::new())]
    outm: String,

    #[arg(long, default_value_t = String::new())]
    outu: String,

    #[arg(long, default_value_t = String::new())]
    outm2: String,

    #[arg(long, default_value_t = String::new())]
    outu2: String,

    #[arg(long)]
    interinput: bool,
}

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

fn main() {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let args = Args::parse();
    duk::run(args);
}
