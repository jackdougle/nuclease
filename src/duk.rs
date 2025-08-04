extern crate bincode;
extern crate needletail;

use crate::kmer_processor::KmerProcessor;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::sync::Arc;
use std::sync::mpsc::channel;
use std::time::Instant;

// TODO:
// - Improve multi-threading
// - Add paired reads support
// - Check out SIMD encoding
// - Check Rust-Bio for all kmer/hashing operations
// - Add memory maximum (argument)
// - Add # of threads to use (argument)

pub fn run(args: crate::Args) {
    let start_time = Instant::now();

    let k = args.k;
    let threshold = args.num_hits;
    let _num_threads = args.threads;
    let ref_path = &args.ref_path;
    let in_path = &args.in_path;
    let outm_path = &args.matched_path;
    let out_path = &args.unmatched_path;
    let bin_kmers_path = &args.bin_kmers_path;
    let _interleaved = args.x;

    let mut kmer_processor = KmerProcessor::new(k, threshold);

    println!("Executing Rust-Duk [{:#?}]\nVersion 1.0.0", args);
    println!(
        "Indexing k-mer table:\t{:.3}",
        start_time.elapsed().as_secs_f32()
    );

    match load_serialized_kmers(bin_kmers_path, &mut kmer_processor) {
        Ok(()) => {
            println!(
                "\nLoaded {} k-mers from {}",
                kmer_processor.ref_kmers.iter().size_hint().0,
                bin_kmers_path
            )
        }
        Err(e) => {
            eprintln!(
                "\nInvalid/no seralized k-mers found: {}\nLoading references from {}",
                e, ref_path
            );

            match get_reference_kmers(ref_path, &mut kmer_processor) {
                Ok(()) => println!(
                    "Added {} from {}",
                    kmer_processor.ref_kmers.iter().size_hint().0,
                    ref_path,
                ),
                Err(e) => {
                    eprintln!("Error loading reference sequences: {}", e);
                    std::process::exit(1);
                }
            };

            match serialize_kmers(bin_kmers_path, &mut kmer_processor) {
                Ok(()) => println!("Saved serialized k-mers to {}", bin_kmers_path),
                Err(e) => eprintln!("Could not serialize reference k-mers: {}", e),
            }
        }
    }

    let indexing_time = start_time.elapsed().as_secs_f32();
    println!("Indexing time:\t\t{:.3} seconds\n", indexing_time);
    println!("Processing reads from {}", in_path);

    match process_reads_crossbeam(&in_path, kmer_processor, &outm_path, &out_path) {
        Ok((mseq_count, mbase_count, useq_count, ubase_count)) => {
            let read_count = mseq_count + useq_count;
            let matched_percent = (mseq_count as f32 / read_count as f32) * 100.0;
            let unmatched_percent = (useq_count as f32 / read_count as f32) * 100.0;

            let base_count = mbase_count + ubase_count;
            let mbase_percent = (mbase_count as f32 / base_count as f32) * 100.0;
            let ubase_percent = (ubase_count as f32 / base_count as f32) * 100.0;

            let end_time = start_time.elapsed().as_secs_f32();
            println!("Processing time:\t{:.3} seconds", end_time - indexing_time);

            println!("\nInput:\t\t\t{} reads", read_count);
            println!(
                "Matches:\t\t{} reads ({:.2}%)\t\t{} bases ({:.2}%)",
                mseq_count, matched_percent, mbase_count, mbase_percent
            );
            println!(
                "Nonmatches:\t\t{} reads ({:.2}%)\t\t{} bases ({:.2}%)",
                useq_count, unmatched_percent, ubase_count, ubase_percent
            );

            println!("\nTime:\t\t\t{:.3} seconds", end_time);

            if read_count >= 10_000 {
                let read_count = read_count as f32 / 1_000 as f32;
                println!(
                    "Reads Processed:\t{:.2}k reads\t\t\t{:.2}k reads/sec",
                    read_count,
                    read_count / end_time
                );
            } else {
                println!(
                    "Reads Processed:\t{} reads\t\t\t{:.2} reads/sec",
                    read_count,
                    read_count as f32 / end_time
                );
            }

            if base_count >= 10_000_000 {
                let base_count = base_count as f32 / 1_000_000.0;
                println!(
                    "Bases Processed:\t{:.2}m bases\t\t\t{:.2}m bases/sec",
                    base_count,
                    base_count / end_time
                );
            } else {
                let base_count = base_count as f32 / 1_000.0;
                println!(
                    "Bases Processed:\t{:.2}k bases\t\t\t{:.2}k bases/sec",
                    base_count,
                    base_count / end_time
                );
            }
        }
        Err(e) => {
            eprintln!("Error processing reads:\n{}", e);
        }
    }
}

fn load_serialized_kmers(
    bin_kmers_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let bin_kmers_file = File::open(bin_kmers_path)?;
    let mut reader = BufReader::new(bin_kmers_file);

    processor.ref_kmers = bincode::decode_from_std_read(&mut reader, bincode::config::standard())?;

    Ok(())
}

fn get_reference_kmers(
    ref_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let mut reader = needletail::parse_fastx_file(ref_path)?;

    while let Some(record) = reader.next() {
        let record = record?;

        processor.process_ref(&record.seq());
    }
    Ok(())
}

fn serialize_kmers(path: &str, processor: &mut KmerProcessor) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    bincode::encode_into_std_write(
        &processor.ref_kmers,
        &mut writer,
        bincode::config::standard(),
    )?;

    println!("Saved binary-encoded k-mers to index file: {}", path);
    Ok(())
}

fn process_reads_crossbeam(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let processor = Arc::new(processor);

    let (tx, rx) = channel(); // for processed outputs
    let chunk_size = 1_000;
    let mut chunk = Vec::with_capacity(chunk_size);

    while let Some(record) = reader.next() {
        let record = record?;
        let id: Arc<str> = Arc::from(str::from_utf8(record.id())?);
        let seq = record.seq().into_owned();
        chunk.push((id, seq));

        if chunk.len() == chunk_size {
            let tx = tx.clone();
            let processor = processor.clone();
            let local_chunk = std::mem::take(&mut chunk);

            rayon::spawn(move || {
                let (matched, unmatched) = process_chunk_crossbeam(local_chunk, &processor);
                tx.send((matched, unmatched)).unwrap();
            });
        }
    }

    if !chunk.is_empty() {
        let tx = tx.clone();
        let processor = processor.clone();
        rayon::spawn(move || {
            let (matched, unmatched) = process_chunk_crossbeam(chunk, &processor);
            tx.send((matched, unmatched)).unwrap();
        });
    }

    drop(tx); // Close channel

    let mut matched_writer = BufWriter::new(File::create(matched_output)?);
    let mut unmatched_writer = BufWriter::new(File::create(unmatched_output)?);
    let mut mseq_count = 0;
    let mut mbase_count = 0;
    let mut useq_count = 0;
    let mut ubase_count = 0;

    for (matched, unmatched) in rx {
        mseq_count += matched.len() as u32;
        for (id, seq) in matched {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            mbase_count += seq.len() as u32;
        }

        useq_count += unmatched.len() as u32;
        for (id, seq) in unmatched {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            ubase_count += seq.len() as u32;
        }
    }

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_chunk_crossbeam(
    chunk: Vec<(Arc<str>, Vec<u8>)>,
    processor: &KmerProcessor,
) -> (Vec<(Arc<str>, Vec<u8>)>, Vec<(Arc<str>, Vec<u8>)>) {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();

    for (id, seq) in chunk {
        if processor.process_read(&seq) {
            matched.push((id, seq));
        } else {
            unmatched.push((id, seq));
        }
    }

    (matched, unmatched)
}
