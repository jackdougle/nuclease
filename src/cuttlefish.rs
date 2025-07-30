extern crate bincode;
extern crate needletail;

use crate::kmer_processor::KmerProcessor;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::sync::Arc;

// TODO:
// - Check out SIMD encoding
// - Check Rust-Bio for all kmer/hashing operations
// - Test BBDuk for speed and memory usage
// - Add Rayon multi-threading

pub fn run(args: crate::Args) {
    let k = args.k;
    let threshold = args.threshold;
    let reference_filename = &args.reference;
    let query_filename = &args.query;
    let matched_filename = &args.matched_path;
    let unmatched_filename = &args.unmatched_path;
    let serialized_filename = &args.serialized_kmers_filename;
    let mut kmer_processor = KmerProcessor::new(k, threshold);

    println!("matched_filename: {}", matched_filename);
    println!("unmatched_filename: {}", unmatched_filename);
    println!(
        "k-mer size: {}\nthreshold: {}\nreferences: {}\nqueries: {}",
        k, threshold, reference_filename, query_filename
    );

    match load_serialized_kmers(serialized_filename, &mut kmer_processor) {
        Ok(()) => {
            println!(
                "Successfully loaded serialized k-mers from path: {}",
                serialized_filename
            )
        }
        Err(e) => {
            eprintln!(
                "Invalid/no seralized k-mers found: {}\nLoading references from path: {}",
                e, reference_filename
            );

            match get_reference_kmers(reference_filename, &mut kmer_processor) {
                Ok(()) => println!(
                    "Loaded reference sequences\nProcessing reads from read path: {}",
                    query_filename
                ),
                Err(e) => {
                    eprintln!("Error loading reference sequences: {}", e);
                    std::process::exit(1);
                }
            };

            match serialize_kmers(serialized_filename, &mut kmer_processor) {
                Ok(()) => println!("Saved serialized k-mers to path: {}", serialized_filename),
                Err(e) => eprintln!("Could not serialize reference k-mers: {}", e),
            }
        }
    }

    // match process_reads(
    //     &query_filename,
    //     kmer_processor,
    //     matched_filename,
    //     unmatched_filename,
    // ) {
    //     Ok(()) => {}
    //     Err(e) => {
    //         eprintln!("Error processing reads: {}", e);
    //         std::process::exit(1);
    //     }
    // }

    match get_read_chunks(
        &query_filename,
        kmer_processor,
        &matched_filename,
        &unmatched_filename,
    ) {
        Ok(()) => {}
        Err(_e) => {}
    }
}

fn load_serialized_kmers(
    bin_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let serialized_kmer_file = File::open(bin_path)?;
    let mut reader = BufReader::new(serialized_kmer_file);

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

fn process_reads(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(), Box<dyn Error>> {
    let matched_file = File::create(matched_output)?;
    let mut matched_writer = BufWriter::new(matched_file);
    let mut matched_count = 0usize;

    let unmatched_file = File::create(unmatched_output)?;
    let mut unmatched_writer = BufWriter::new(unmatched_file);
    let mut unmatched_count = 0usize;

    let mut reader = needletail::parse_fastx_file(reads_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let id = str::from_utf8(record.id())?;
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

// fn get_read_chunks(
//     reads_path: &str,
//     processor: KmerProcessor,
//     matched_output: &str,
//     unmatched_output: &str,
// ) -> Result<(), Box<dyn Error>> {
//     let mut reader = needletail::parse_fastx_file(reads_path)?;
//     let mut chunk = Vec::new();
//     let safe_processor = Arc::from(processor);

//     let matched_file = File::create(matched_output)?;
//     let mut matched_writer = BufWriter::new(matched_file);
//     let mut matched_count: u32 = 0;

//     let unmatched_file = File::create(unmatched_output)?;
//     let mut unmatched_writer = BufWriter::new(unmatched_file);
//     let mut unmatched_count: u32 = 0;

//     while let Some(record) = reader.next() {
//         let record = record?;
//         let id = String::from_utf8(Vec::from(record.id()))?;
//         let seq = record.seq().into_owned();

//         chunk.push((id, seq));

//         if chunk.len() == 10_000 {
//             match process_reads_parellel(
//                 chunk,
//                 safe_processor.clone(),
//                 &mut matched_writer,
//                 &mut unmatched_writer,
//                 &mut matched_count,
//                 &mut unmatched_count,
//             ) {
//                 Ok(()) => {}
//                 Err(e) => {
//                     eprintln!("Processing reads unsuccessful: {}", e);
//                     std::process::exit(1);
//                 }
//             }
//             chunk = Vec::new();
//         }
//     }

//     if chunk.len() > 0 && chunk.len() < 10_000 {
//         match process_reads_parellel(
//             chunk,
//             safe_processor.clone(),
//             &mut matched_writer,
//             &mut unmatched_writer,
//             &mut matched_count,
//             &mut unmatched_count,
//         ) {
//             Ok(()) => {}
//             Err(e) => {
//                 eprintln!("Processing reads unsuccessful: {}", e);
//                 std::process::exit(1);
//             }
//         }
//     }

//     println!("Wrote {} matched reads to out/matched.fa!", matched_count);
//     println!(
//         "Wrote {} unmatched reads to out/unmatched.fa!",
//         unmatched_count
//     );

//     Ok(())
// }

// fn process_reads_parellel(
//     reads_chunk: Vec<(String, Vec<u8>)>,
//     processor: Arc<KmerProcessor>,
//     matched_writer: &mut BufWriter<File>,
//     unmatched_writer: &mut BufWriter<File>,
//     matched_count: &mut u32,
//     unmatched_count: &mut u32,
// ) -> Result<(), Box<dyn Error>> {
//     for (id, seq) in reads_chunk.iter() {
//         let is_match = processor.process_read(seq.as_ref());

//         if is_match {
//             writeln!(matched_writer, ">{}", id)?;
//             matched_writer.write_all(seq)?;
//             writeln!(matched_writer)?;
//             *matched_count += 1;
//         } else {
//             writeln!(unmatched_writer, ">{}", id)?;
//             unmatched_writer.write_all(seq)?;
//             writeln!(unmatched_writer)?;
//             *unmatched_count += 1;
//         }
//     }

//     Ok(())
// }

use std::sync::Mutex;

fn get_read_chunks(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let mut chunk = Vec::new();

    let processor = Arc::new(processor);

    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_output)?)));
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_output)?)));

    let matched_count = Arc::new(Mutex::new(0u32));
    let unmatched_count = Arc::new(Mutex::new(0u32));

    let mut chunks: Vec<Vec<(String, Vec<u8>)>> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8(Vec::from(record.id()))?;
        let seq = record.seq().into_owned();

        chunk.push((id, seq));

        if chunk.len() == 10_000 {
            chunks.push(chunk);
            chunk = Vec::new();
        }
    }

    if !chunk.is_empty() {
        chunks.push(chunk);
    }

    chunks.into_par_iter().try_for_each(|chunk| {
        process_reads_parallel(
            chunk,
            processor.clone(),
            matched_writer.clone(),
            unmatched_writer.clone(),
            matched_count.clone(),
            unmatched_count.clone(),
        )
    })?;

    println!(
        "Wrote {} matched reads to {}!",
        matched_count.lock().unwrap(),
        matched_output
    );
    println!(
        "Wrote {} unmatched reads to {}!",
        unmatched_count.lock().unwrap(),
        unmatched_output
    );

    Ok(())
}

fn process_reads_parallel(
    reads_chunk: Vec<(String, Vec<u8>)>,
    processor: Arc<KmerProcessor>,
    matched_writer: Arc<Mutex<BufWriter<File>>>,
    unmatched_writer: Arc<Mutex<BufWriter<File>>>,
    matched_count: Arc<Mutex<u32>>,
    unmatched_count: Arc<Mutex<u32>>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let results: Vec<(bool, String, Vec<u8>)> = reads_chunk
        .into_par_iter()
        .map(|(id, seq)| {
            let is_match = processor.process_read(&seq);
            (is_match, id, seq)
        })
        .collect();

    let mut matched_writer = matched_writer.lock().unwrap();
    let mut unmatched_writer = unmatched_writer.lock().unwrap();
    let mut matched_count = matched_count.lock().unwrap();
    let mut unmatched_count = unmatched_count.lock().unwrap();

    for (is_match, id, seq) in results {
        if is_match {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            *matched_count += 1;
        } else {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            *unmatched_count += 1;
        }
    }

    Ok(())
}
