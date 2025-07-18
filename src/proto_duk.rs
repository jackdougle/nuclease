extern crate bincode;
extern crate indicatif;
extern crate needletail;
extern crate rayon;

use crate::kmer::canonical_kmer;

// TODO:
// - Add tests for program speed and memory usage
// - Replace slow logic using Strings and Vecs with faster logic using Bytes and BytesMut
// - Migrate repeated k-mer logic to more efficient logic in kmers.rs
// - Replace Vecs with Arc<[T]> for thread-safe access to Vec<T>
// - Add a progress bar
// - Add a way to specify the number of threads to use
// - Add a way to specify the number of reads to process
// - Add a way to specify the number of k-mers to process
// - Add a way to specify the number of reference sequences to process
// - Add a way to specify the number of query sequences to process

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

pub fn run(args: crate::Args) {
    let k = args.k;
    let threshold = args.threshold;
    let canonical = args.canonical;
    let reference_filename = &args.reference;
    let query_filename = &args.query;
    let matched_filename = &args.matched_path;
    println!("matched_filename: {}", matched_filename);
    let unmatched_filename = &args.unmatched_path;
    println!("unmatched_filename: {}", unmatched_filename);
    let serialized_kmers_filename = &args.serialized_kmers_filename;

    println!(
        "k-mer size: {}\nthreshold: {}\ncanonical: {}\nreferences: {}\nqueries: {}",
        k, threshold, canonical, reference_filename, query_filename
    );

    // Try to load pre-built k-mer index, or build it if it doesn't exist
    let ref_kmers = match load_kmer_index(serialized_kmers_filename, k) {
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
            let ref_seqs = match load_reference_sequences(reference_filename) {
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
            if let Err(e) = save_kmer_index(&kmers, serialized_kmers_filename) {
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
            let _ =
                write_results_with_ids(&matched, &unmatched, matched_filename, unmatched_filename);
        }
        Err(e) => {
            eprintln!("Error processing reads: {}", e);
            std::process::exit(1);
        }
    }
}

/// Loads k-mer index from disk for fast startup
fn load_kmer_index(path: &str, k: usize) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::BufReader;

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // Deserialize k-mer set from binary format
    let kmers: HashSet<String> =
        bincode::decode_from_std_read(&mut reader, bincode::config::standard())?;
    if kmers.iter().next().unwrap().len() != k {
        println!(
            "k-mer size mismatch: {} != {}, reference binary file invalid",
            kmers.iter().next().unwrap().len(),
            k
        );
        std::process::exit(1);
    }

    Ok(kmers)
}

fn get_reference_kmers(
    ref_seqs: &HashMap<Arc<str>, Arc<str>>,
    k: usize,
    canonical: bool,
) -> HashSet<String> {
    let mut ref_kmers: HashSet<String> = HashSet::new();

    for (_name, seq) in ref_seqs {
        for x in 0..=seq.len() - k {
            let temp = seq[x..x + k].to_string();

            if canonical {
                ref_kmers.insert(canonical_kmer(&temp));
            } else {
                ref_kmers.insert(temp);
            }
        }
    }

    ref_kmers
}

fn load_reference_sequences(
    path: &str,
) -> Result<HashMap<Arc<str>, Arc<str>>, Box<dyn std::error::Error>> {
    let mut ref_seqs: HashMap<Arc<str>, Arc<str>> = HashMap::new();

    // needletail automatically detects file format (FASTA/FASTQ) and handles compression
    let mut reader = needletail::parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;

        // Convert byte slices to strings (needletail provides data as &[u8])
        let id = std::str::from_utf8(record.id()).unwrap();
        let seq = String::from_utf8(record.seq().to_vec()).unwrap();

        ref_seqs.insert(Arc::from(id), Arc::from(seq.as_str()));
    }

    Ok(ref_seqs)
}

/// Writes results to output files
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
                read_kmers.insert(canonical_kmer(kmer));
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
fn save_kmer_index(kmers: &HashSet<String>, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::BufWriter;

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Serialize k-mer set to binary format
    bincode::encode_into_std_write(kmers, &mut writer, bincode::config::standard())?;

    println!("Saved {} k-mers to index file: {}", kmers.len(), path);
    Ok(())
}
