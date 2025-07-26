extern crate bincode;
extern crate needletail;

use crate::kmer_processor::KmerProcessor;
use std::fs::File;
use std::io::{BufWriter, Write};

// TODO:
// - Test BBDuk for speed and memory usage
// - Add Rayon multi-threading
// - Check Rust-Bio for k-mer utilities
// - Add error handling for k-mer index loading

pub fn run(args: crate::Args) {
    let k = args.k;
    let threshold = args.threshold;
    let reference_filename = &args.reference;
    let query_filename = &args.query;
    let matched_filename = &args.matched_path;
    let unmatched_filename = &args.unmatched_path;

    let mut kmer_processor = KmerProcessor::new(k, threshold);

    println!("matched_filename: {}", matched_filename);
    println!("unmatched_filename: {}", unmatched_filename);
    println!(
        "k-mer size: {}\nthreshold: {}\nreferences: {}\nqueries: {}",
        k, threshold, reference_filename, query_filename
    );

    match load_reference_sequences(reference_filename, &mut kmer_processor) {
        Ok(_ok) => {
            println!(
                "Loaded reference sequences\nProcessing reads from read path: {}",
                query_filename
            );
        }
        Err(e) => {
            eprintln!("Error loading reference sequences: {}", e);
            std::process::exit(1);
        }
    };

    match process_reads_and_write_streaming(
        &query_filename,
        kmer_processor,
        matched_filename,
        unmatched_filename,
    ) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error processing reads: {}", e);
            std::process::exit(1);
        }
    }
}

fn load_reference_sequences(
    path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = needletail::parse_fastx_file(path)?;

    while let Some(record) = reader.next() {
        let record = record?;

        processor.process_ref(&record.seq());
    }
    Ok(())
}

fn process_reads_and_write_streaming(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let matched_file = File::create(matched_output)?;
    let mut matched_writer = BufWriter::new(matched_file);
    let mut matched_count = 0usize;

    let unmatched_file = File::create(unmatched_output)?;
    let mut unmatched_writer = BufWriter::new(unmatched_file);
    let mut unmatched_count = 0usize;

    let mut reader = needletail::parse_fastx_file(reads_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8_lossy(record.id());
        let seq = &record.seq();

        let is_match = processor.process_read(seq);

        if is_match {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(seq)?;
            writeln!(matched_writer)?;
            matched_count += 1;
        } else {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(seq)?;
            writeln!(unmatched_writer)?;
            unmatched_count += 1;
        }
    }

    println!("Wrote {} matched reads to out/matched.fa!", matched_count);
    println!(
        "Wrote {} unmatched reads to out/unmatched.fa!",
        unmatched_count
    );

    Ok(())
}
