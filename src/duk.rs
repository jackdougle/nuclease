use crate::kmer_processor::KmerProcessor;
use bincode::{config, decode_from_std_read, encode_into_std_write};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::mem::take;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::thread;
use std::time::Instant;

// TODO:
// - Add paired reads support
// - Add k-mer length metadata to serialized file and error checking
// - Add memory maximum (argument)
// - Add # of threads to use (argument)
// - Do documentation (vvv)
// - Github repo: follow crates.io syntax
// - Try to upload crate to crates.io and docs.rs
// - Consider pull request to Rust-Bio (still going?)

pub fn run(args: crate::Args) {
    let start_time = Instant::now();

    let k = args.k;
    let min_hits = args.minhits;
    let _num_threads = args.threads;
    let _max_memory = args.maxmem;
    let ref_path = args.r#ref;
    let in_path = args.r#in;
    let in2_path = args.in2;
    let outm_path = args.outm;
    let outu_path = args.outu;
    let outm2_path = args.outm2;
    let outu2_path = args.outu2;
    let bin_kmers_path = &args.binref;
    let interleaved_input = args.interinput;

    let mut kmer_processor = KmerProcessor::new(k, min_hits);

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

            match get_reference_kmers(&ref_path, &mut kmer_processor) {
                Ok(()) => println!(
                    "Added {} from {}",
                    kmer_processor.ref_kmers.iter().size_hint().0,
                    ref_path,
                ),
                Err(e) => {
                    eprintln!("\nError loading reference sequences: {}", e);
                    std::process::exit(1);
                }
            };

            match serialize_kmers(bin_kmers_path, &mut kmer_processor) {
                Ok(()) => println!("Saved serialized k-mers to {}", bin_kmers_path),
                Err(e) => eprintln!("\nCould not serialize reference k-mers: {}", e),
            }
        }
    }

    let indexing_time = start_time.elapsed().as_secs_f32();
    println!("Indexing time:\t\t{:.3} seconds\n", indexing_time);
    println!("Processing reads from {}", in_path);

    match process_reads(
        in_path,
        in2_path,
        kmer_processor,
        &outm_path,
        &outu_path,
        &outm2_path,
        &outu2_path,
        interleaved_input,
    ) {
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
                "Matches:\t\t{} reads ({:.2}%)   \t\t{} bases ({:.2}%)",
                mseq_count, matched_percent, mbase_count, mbase_percent
            );
            println!(
                "Nonmatches:\t\t{} reads ({:.2}%)\t\t{} bases ({:.2}%)\n",
                useq_count, unmatched_percent, ubase_count, ubase_percent
            );

            println!(
                "Time:\t\t\t{:.3} seconds",
                start_time.elapsed().as_secs_f32()
            );

            if read_count >= 1_000_000 {
                let read_count = read_count as f32 / 1_000_000.0;
                println!(
                    "Reads Processed:\t{:.2}m reads\t\t\t{:.2}m reads/sec",
                    read_count,
                    read_count / end_time
                );
            } else if read_count >= 10_000 {
                let read_count = read_count as f32 / 1_000.0;
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

            if base_count >= 1_000_000 {
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
            // TODO: Add error handling for missing file inputs (needs both >0 matched and unmatched files right now)
            eprintln!("\nError processing reads:\n{}", e);
        }
    }
}

fn load_serialized_kmers(
    bin_kmers_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let bin_kmers_file = File::open(bin_kmers_path)?;
    let mut reader = BufReader::new(bin_kmers_file);

    processor.ref_kmers = decode_from_std_read(&mut reader, config::standard())?;

    Ok(())
}

fn get_reference_kmers(
    ref_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let mut reader = parse_fastx_file(ref_path)?;

    while let Some(record) = reader.next() {
        let record = record?;

        processor.process_ref(&record.seq());
    }
    Ok(())
}

fn serialize_kmers(path: &str, processor: &mut KmerProcessor) -> Result<(), Box<dyn Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    encode_into_std_write(&processor.ref_kmers, &mut writer, config::standard())?;

    println!("Saved binary-encoded k-mers to index file: {}", path);
    Ok(())
}

#[derive(PartialEq, Clone, Copy, Debug, Default)]
enum ProcessMode {
    #[default]
    Unpaired,
    Paired,
    PairedInInterOut,
    InterInPairedOut,
    Interleaved,
}

fn detect_mode(
    reads2_path: &str,
    matched2_path: &str,
    unmatched2_path: &str,
    interleaved_input: bool,
) -> ProcessMode {
    if !reads2_path.is_empty() {
        if matched2_path.is_empty() && unmatched2_path.is_empty() {
            println!(
                "Forcing interleaved output because paired input was specified for single output files"
            );
            ProcessMode::PairedInInterOut
        } else {
            assert!(
                !matched2_path.is_empty(),
                "Please add a second matched output path using: --outm2 <file>"
            );
            assert!(
                !unmatched2_path.is_empty(),
                "Please add a second unmatched output path using: --outu2 <file>"
            );
            println!("Input and output is processed as paired");
            ProcessMode::Paired
        }
    } else if interleaved_input {
        if matched2_path.is_empty() && unmatched2_path.is_empty() {
            println!("Input and output is processed as interleaved");
            ProcessMode::Interleaved
        } else {
            assert!(
                !matched2_path.is_empty(),
                "Please add a second matched output path using: --outm2 <file>"
            );
            assert!(
                !unmatched2_path.is_empty(),
                "Please add a second unmatched output path using: --outu2 <file>"
            );
            println!("Processing interleaved input and paired output");
            ProcessMode::InterInPairedOut
        }
    } else if !matched2_path.is_empty() && !unmatched2_path.is_empty() {
        println!(
            "Forcing interleaved input because paired output was specified for a single input file"
        );
        ProcessMode::InterInPairedOut
    } else {
        println!("Input and output are processed as unpaired");
        ProcessMode::Unpaired
    }
}

fn process_reads(
    reads_path: String,
    reads2_path: String,
    processor: KmerProcessor,
    matched_path: &str,
    unmatched_path: &str,
    matched2_path: &str,
    unmatched2_path: &str,
    interleaved_input: bool,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let processor = Arc::new(processor);
    let process_mode = detect_mode(
        &reads2_path,
        matched2_path,
        unmatched2_path,
        interleaved_input,
    );

    let (sender, receiver): (
        Sender<Vec<(bool, Vec<u8>, Vec<u8>)>>,
        Receiver<Vec<(bool, Vec<u8>, Vec<u8>)>>,
    ) = channel();

    // Reader thread
    let sender2 = sender.clone();
    let reader_thread = thread::spawn(move || -> Result<(), Box<dyn Error + Send + Sync>> {
        let mut reader = parse_fastx_file(reads_path)?;

        const CHUNK_SIZE: usize = 2_000;
        let mut chunk = Vec::with_capacity(CHUNK_SIZE);

        if process_mode == ProcessMode::Unpaired
            || process_mode == ProcessMode::Interleaved
            || process_mode == ProcessMode::InterInPairedOut
        {
            while let Some(record) = reader.next() {
                let record = record?;
                chunk.push((record.id().to_vec(), record.seq().to_vec()));

                if chunk.len() == CHUNK_SIZE {
                    let processor = processor.clone();
                    let sender = sender2.clone();
                    let local_chunk = take(&mut chunk);

                    // Process chunk in parallel
                    rayon::spawn(move || {
                        let processed: Vec<(bool, Vec<u8>, Vec<u8>)> = local_chunk
                            .into_par_iter()
                            .map(|(id, seq)| {
                                let has_match = processor.process_read(&seq);
                                (has_match, id, seq)
                            })
                            .collect();

                        let _ = sender.send(processed);
                    });
                }
            }

            if !chunk.is_empty() {
                let processed: Vec<(bool, Vec<u8>, Vec<u8>)> = chunk
                    .into_par_iter()
                    .map(|(id, seq)| {
                        let has_match = processor.process_read(&seq);
                        (has_match, id, seq)
                    })
                    .collect();

                let _ = sender2.send(processed);
            }
        }

        if process_mode == ProcessMode::PairedInInterOut || process_mode == ProcessMode::Paired {
            let mut reader2 = parse_fastx_file(reads2_path)?;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            while let (Some(record1), Some(record2)) = (reader.next(), reader2.next()) {
                let record1 = record1?;
                let record2 = record2?;

                chunk.push((record1.id().to_vec(), record1.seq().to_vec()));
                chunk.push((record2.id().to_vec(), record2.seq().to_vec()));

                if chunk.len() == CHUNK_SIZE {
                    let processor = processor.clone();
                    let sender = sender2.clone();
                    let local_chunk = take(&mut chunk);

                    rayon::spawn(move || {
                        let processed: Vec<(bool, Vec<u8>, Vec<u8>)> = local_chunk
                            .into_par_iter()
                            .map(|(id, seq)| {
                                let has_match = processor.process_read(&seq);
                                // decide logic: both must match? either?
                                (has_match, id, seq)
                            })
                            .collect();

                        let _ = sender.send(processed);
                    });
                }
            }

            if !chunk.is_empty() {
                let processed: Vec<(bool, Vec<u8>, Vec<u8>)> = chunk
                    .into_par_iter()
                    .map(|(id, seq)| {
                        let has_match = processor.process_read(&seq);
                        (has_match, id, seq)
                    })
                    .collect();

                let _ = sender2.send(processed);
            }
        }

        Ok(())
    });

    drop(sender);

    let mut matched_writer = BufWriter::new(File::create(matched_path)?);
    let mut unmatched_writer = BufWriter::with_capacity(4_000_000, File::create(unmatched_path)?);

    let matched_count = Arc::new(AtomicU32::new(0));
    let matched_bases = Arc::new(AtomicU32::new(0));
    let unmatched_count = Arc::new(AtomicU32::new(0));
    let unmatched_bases = Arc::new(AtomicU32::new(0));

    if process_mode == ProcessMode::Unpaired {
        for batch in receiver {
            for (has_match, id, seq) in batch {
                if has_match {
                    matched_writer.write_all(b">")?;
                    matched_writer.write_all(&id)?;
                    matched_writer.write_all(b"\n")?;
                    matched_writer.write_all(&seq)?;
                    matched_writer.write_all(b"\n")?;

                    matched_count.fetch_add(1, Ordering::Relaxed);
                    matched_bases.fetch_add(seq.len() as u32, Ordering::Relaxed);
                } else {
                    unmatched_writer.write_all(b">")?;
                    unmatched_writer.write_all(&id)?;
                    unmatched_writer.write_all(b"\n")?;
                    unmatched_writer.write_all(&seq)?;
                    unmatched_writer.write_all(b"\n")?;

                    unmatched_count.fetch_add(1, Ordering::Relaxed);
                    unmatched_bases.fetch_add(seq.len() as u32, Ordering::Relaxed);
                }
            }
        }
    } else if process_mode == ProcessMode::Interleaved
        || process_mode == ProcessMode::PairedInInterOut
    {
        for batch in receiver {
            for i in (0..batch.len() - 1).step_by(2) {
                let has_match = batch[i].0 || batch[i + 1].0;
                if has_match {
                    matched_writer.write_all(b">")?;
                    matched_writer.write_all(&batch[i].1)?;
                    matched_writer.write_all(b"\n")?;
                    matched_writer.write_all(&batch[i].2)?;
                    matched_writer.write_all(b"\n")?;

                    matched_writer.write_all(b">")?;
                    matched_writer.write_all(&batch[i + 1].1)?;
                    matched_writer.write_all(b"RC\n")?;
                    matched_writer.write_all(&batch[i + 1].2)?;
                    matched_writer.write_all(b"\n")?;

                    matched_count.fetch_add(1, Ordering::Relaxed);
                    matched_bases.fetch_add(batch[i].2.len() as u32, Ordering::Relaxed);
                } else {
                    unmatched_writer.write_all(b">")?;
                    unmatched_writer.write_all(&batch[i].1)?;
                    unmatched_writer.write_all(b"\n")?;
                    unmatched_writer.write_all(&batch[i].2)?;
                    unmatched_writer.write_all(b"\n")?;

                    unmatched_writer.write_all(b">")?;
                    unmatched_writer.write_all(&batch[i + 1].1)?;
                    unmatched_writer.write_all(b"RC\n")?;
                    unmatched_writer.write_all(&batch[i + 1].2)?;
                    unmatched_writer.write_all(b"\n")?;

                    unmatched_count.fetch_add(1, Ordering::Relaxed);
                    unmatched_bases.fetch_add(batch[i].2.len() as u32, Ordering::Relaxed);
                }
            }
        }
    } else if process_mode == ProcessMode::InterInPairedOut || process_mode == ProcessMode::Paired {
        let mut matched2_writer: BufWriter<File> = BufWriter::new(File::create(matched2_path)?);
        let mut unmatched2_writer: BufWriter<File> = BufWriter::new(File::create(unmatched2_path)?);

        for batch in receiver {
            for i in (0..batch.len() - 1).step_by(2) {
                let has_match = batch[i].0 || batch[i + 1].0;
                if has_match {
                    matched_writer.write_all(b">")?;
                    matched_writer.write_all(&batch[i].1)?;
                    matched_writer.write_all(b"\n")?;
                    matched_writer.write_all(&batch[i].2)?;
                    matched_writer.write_all(b"\n")?;

                    matched2_writer.write_all(b">")?;
                    matched2_writer.write_all(&batch[i + 1].1)?;
                    matched2_writer.write_all(b"\n")?;
                    matched2_writer.write_all(&batch[i + 1].2)?;
                    matched2_writer.write_all(b"\n")?;

                    matched_count.fetch_add(1, Ordering::Relaxed);
                    matched_bases.fetch_add(batch[i].2.len() as u32, Ordering::Relaxed);
                } else {
                    unmatched_writer.write_all(b">")?;
                    unmatched_writer.write_all(&batch[i].1)?;
                    unmatched_writer.write_all(b"\n")?;
                    unmatched_writer.write_all(&batch[i].2)?;
                    unmatched_writer.write_all(b"\n")?;

                    unmatched2_writer.write_all(b">")?;
                    unmatched2_writer.write_all(&batch[i + 1].1)?;
                    unmatched2_writer.write_all(b"\n")?;
                    unmatched2_writer.write_all(&batch[i + 1].2)?;
                    unmatched2_writer.write_all(b"\n")?;

                    unmatched_count.fetch_add(1, Ordering::Relaxed);
                    unmatched_bases.fetch_add(batch[i].2.len() as u32, Ordering::Relaxed);
                }
            }
        }

        matched2_writer.flush()?;
        unmatched2_writer.flush()?;
    }

    matched_writer.flush()?;
    unmatched_writer.flush()?;

    reader_thread.join().unwrap()?;

    Ok((
        matched_count.load(Ordering::Relaxed),
        matched_bases.load(Ordering::Relaxed),
        unmatched_count.load(Ordering::Relaxed),
        unmatched_bases.load(Ordering::Relaxed),
    ))
}
