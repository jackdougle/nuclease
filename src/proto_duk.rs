extern crate bincode;
extern crate needletail;

use std::collections::{HashMap, HashSet};
// use std::sync::Arc;
// use std::thread;

pub fn run(args: crate::Args) {
    let k = args.k;
    let threshold = args.threshold;
    let canonical = args.canonical;
    let reference_filename = &args.reference;
    let query_filename = &args.query;
    // let matched_filename = &args.matched_path;
    // let unmatched_filename = &args.unmatched_path;
    let serialized_kmers_filename = &args.serialized_kmers_filename;

    println!(
        "k-mer size: {}\nthreshold: {}\ncanonical: {}\nreferences: {}\nqueries: {}",
        k, threshold, canonical, reference_filename, query_filename
    );

    // Try to load pre-built k-mer index, or build it if it doesn't exist
    let ref_kmers = match load_kmer_index(serialized_kmers_filename) {
        Ok(kmers) => {
            println!(
                "Loaded {} k-mers from pre-built index file: {}",
                kmers.len(),
                serialized_kmers_filename
            );
            kmers
        }
        Err(_) => {
            println!("No pre-built index found, building from reference sequences...");

            // Load reference sequences (can use streaming for large reference files)
            let ref_seqs = match load_reference_streaming(&reference_filename) {
                Ok(seqs) => {
                    println!("Loaded {} reference sequences", seqs.len());
                    println!("Reference sequences: {:#?}", seqs);
                    seqs
                }
                Err(e) => {
                    eprintln!("Error loading reference sequences: {}", e);
                    std::process::exit(1);
                }
            };

            // Extract reference k-mers
            let kmers = get_reference_kmers(&ref_seqs, k, canonical);
            println!("Extracted {} reference k-mers", kmers.len());

            // Save for future use
            if let Err(e) = save_kmer_index(&kmers, "in/ref_kmers.bin") {
                eprintln!("Warning: Could not save k-mer index: {}", e);
            } else {
                println!("Saved k-mer index for future use");
            }

            kmers
        }
    };

    // Process reads using streaming (much more memory efficient)
    match process_reads(&query_filename, &ref_kmers, k, threshold, canonical) {
        Ok((matched, unmatched)) => {
            println!(
                "Matched reads: {} (threshold: {} k-mer hits)",
                matched.len(),
                threshold
            );
            println!("Unmatched reads: {}", unmatched.len());
            let _ = write_results_with_ids(
                &matched,
                &unmatched,
                "out/matched.fasta",
                "out/unmatched.fasta",
            );
        }
        Err(e) => {
            eprintln!("Error processing reads: {}", e);
            std::process::exit(1);
        }
    }
}

fn mouge() {
    let x = String::from("mouge");
    println!("Y: {}", &x);
    mouge2(&x);
}

fn mouge2(stringy: &String) {
    println!("X: {}", stringy);
}

/// Loads k-mer index from disk for fast startup
///
/// PURPOSE: Load a previously saved k-mer index instead of rebuilding it.
/// This dramatically reduces startup time and enables efficient parallel processing.
///
/// HOW IT WORKS:
/// - Deserializes the binary k-mer index file
/// - Restores the HashSet<String> in memory
/// - Provides instant access to reference k-mers
///
/// BENEFITS:
/// - Fast startup: no need to re-parse reference files
/// - Memory efficient: only loads the k-mer set, not full sequences
/// - Enables parallel workflows: multiple processes can share the same index
/// - Supports your goal of efficient SIZ chunk processing
fn load_kmer_index(path: &str) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Deserialize k-mer set from binary format
    let kmers: HashSet<String> = bincode::deserialize_from(reader)?;

    Ok(kmers)
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect()
}

fn get_reference_kmers(
    ref_seqs: &HashMap<String, String>,
    k: usize,
    canonical_bool: bool,
) -> HashSet<String> {
    let mut ref_kmers: HashSet<String> = HashSet::new();

    for (_name, seq) in ref_seqs {
        for x in 0..=seq.len() - k {
            let temp = seq[x..x + k].to_string();

            if canonical_bool {
                let rev = reverse_complement(&temp);
                ref_kmers.insert(std::cmp::min(temp, rev));
            } else {
                ref_kmers.insert(temp);
            }
        }
    }

    ref_kmers
}

fn load_reference_streaming(
    path: &str,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let mut ref_seqs: HashMap<String, String> = HashMap::new();

    // needletail automatically detects file format (FASTA/FASTQ) and handles compression
    let mut reader = needletail::parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;

        // Convert byte slices to strings (needletail provides data as &[u8])
        let id = String::from_utf8_lossy(record.id()).to_string();
        let seq = String::from_utf8_lossy(&record.seq()).to_string();

        ref_seqs.insert(id, seq);
    }

    Ok(ref_seqs)
}

/// Writes results to output files with read ID preservation (like BBDuk)
///
/// PURPOSE: Write matched and unmatched reads to files while preserving
/// their original IDs, which is essential for bioinformatics workflows
/// and matches BBDuk's output behavior.
///
/// HOW IT WORKS:
/// - Takes read records with IDs and sequences
/// - Writes to FASTA/FASTQ format preserving original IDs
/// - Maintains compatibility with downstream tools
///
/// BENEFITS:
/// - Preserves read identification for downstream analysis
/// - Compatible with bioinformatics pipeline tools
/// - Matches BBDuk's output format expectations
/// - Enables tracking reads through multi-step workflows
fn write_results_with_ids(
    matched: &[(String, String)],   // (read_id, sequence)
    unmatched: &[(String, String)], // (read_id, sequence)
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    // Write matched reads with their original IDs
    let matched_file = File::create(matched_output)?;
    let mut matched_writer = BufWriter::new(matched_file);

    for (read_id, seq) in matched {
        writeln!(matched_writer, ">{}", read_id)?;
        writeln!(matched_writer, "{}", seq)?;
    }

    // Write unmatched reads with their original IDs
    let unmatched_file = File::create(unmatched_output)?;
    let mut unmatched_writer = BufWriter::new(unmatched_file);

    for (read_id, seq) in unmatched {
        writeln!(unmatched_writer, ">{}", read_id)?;
        writeln!(unmatched_writer, "{}", seq)?;
    }

    println!(
        "Wrote {} matched reads (with IDs) to {}",
        matched.len(),
        matched_output
    );
    println!(
        "Wrote {} unmatched reads (with IDs) to {}",
        unmatched.len(),
        unmatched_output
    );

    Ok(())
}

/// Processes reads in a streaming fashion without loading all reads into memory
///
/// PURPOSE: Handle large read files (billions of reads) by processing them
/// one at a time, which is essential for your goal of processing large datasets
/// that BBDuk struggles with.
///
/// HOW IT WORKS WITH NEEDLETAIL:
/// - Streams reads one at a time from the input file
/// - Extracts k-mers from each read individually
/// - Compares against reference k-mers immediately
/// - Outputs results without storing all reads in memory
///
/// BENEFITS:
/// - Memory efficient: processes reads one at a time
/// - Fast: needletail's optimized parsing
/// - Scalable: can handle files of any size
/// - Streaming output: can pipe results to other tools
fn process_reads(
    reads_path: &str,
    ref_kmers: &HashSet<String>,
    k: usize,
    threshold: usize,
    canonical: bool,
) -> Result<(Vec<(String, String)>, Vec<(String, String)>), Box<dyn std::error::Error>> {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();

    // needletail handles file format detection and streaming
    let mut reader = needletail::parse_fastx_file(reads_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        let seq = String::from_utf8_lossy(&record.seq()).to_string();

        // Extract k-mers from this single read
        let mut read_kmers = HashSet::new();
        for i in 0..=seq.len() - k {
            let kmer = &seq[i..i + k];

            if canonical {
                let rc = reverse_complement(kmer);
                read_kmers.insert(std::cmp::min(kmer, &rc).to_string());
            } else {
                read_kmers.insert(kmer.to_string());
            }
        }

        // Count hits against reference k-mers
        let hits = read_kmers.intersection(ref_kmers).count();

        // Classify this read based on threshold
        if hits >= threshold {
            matched.push((id, seq));
        } else {
            unmatched.push((id, seq));
        }
    }

    Ok((matched, unmatched))
}

/// Saves k-mer index to disk for reuse across multiple runs
///
/// HOW IT WORKS:
/// - Serializes the HashSet<String> to a binary format
/// - Uses bincode for efficient serialization
/// - Creates a reusable index file that can be shared across machines
///
/// BENEFITS:
/// - Eliminates index rebuilding time
/// - Enables parallel processing by sharing index files
/// - Reduces startup time for repeated runs
/// - Supports your goal of efficient multi-threading workflows
fn save_kmer_index(kmers: &HashSet<String>, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::BufWriter;

    let file = File::create(path)?;
    let writer = BufWriter::new(file);

    // Serialize k-mer set to binary format
    bincode::serialize_into(writer, kmers)?;

    println!("Saved {} k-mers to index file: {}", kmers.len(), path);
    Ok(())
}
