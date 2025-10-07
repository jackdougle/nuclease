use crate::kmer_ops::KmerProcessor;
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
use std::time::Instant;
use std::{fs, io, u32};
use std::{thread, usize};

pub fn run(args: crate::Args, start_time: Instant) -> io::Result<()> {
    let available_threads = num_cpus::get();

    // First set of variables: reference indexing parameters
    let num_threads = args
        .threads
        .unwrap_or(available_threads)
        .min(available_threads);
    let k = args.k.unwrap_or(21);
    let min_hits = args.minhits.unwrap_or(1);
    let ordered_output = args.order;

    let ref_path = args.r#ref;
    let bin_kmers_path = &args.binref.unwrap_or_default();
    let new_bin_kmers_path = &args.saveref.unwrap_or_default();

    let mut kmer_processor = KmerProcessor::new(k, min_hits);

    match load_serialized_kmers(bin_kmers_path, &mut kmer_processor) {
        Ok(()) => {
            println!(
                "\nLoaded {} k-mers from {}",
                kmer_processor.ref_kmers.iter().size_hint().0 - 1,
                bin_kmers_path
            )
        }
        Err(e) => {
            eprintln!("\nInvalid serialized reference file: {}", e);
            println!("Loading ref k-mers from {}", ref_path);

            match get_reference_kmers(&ref_path, &mut kmer_processor) {
                Ok(()) => println!(
                    "Added {} from {}",
                    kmer_processor.ref_kmers.iter().size_hint().0 - 1,
                    ref_path,
                ),
                Err(e) => {
                    eprintln!("\nError loading reference sequences: {}", e);
                    std::process::exit(1);
                }
            };

            match serialize_kmers(&new_bin_kmers_path, &mut kmer_processor) {
                Ok(()) => println!("Saved serialized k-mers to {}", new_bin_kmers_path),
                Err(e) => eprintln!("\nCould not serialize reference k-mers: {}", e),
            }
        }
    }

    let indexing_time = start_time.elapsed().as_secs_f32();
    println!("Indexing time:\t\t{:.3} seconds\n", indexing_time);

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Could not build Rayon Pool with specified thread amount");

    let process_mode = detect_mode(&args.in2, &args.outm2, &args.outu2, args.interinput);

    // Second set of variables: read processing I/O
    let in_path = args.r#in;
    let in2_path = args.in2.unwrap_or_default();

    let outm_path = args.outm.unwrap_or(String::from("/dev/null"));
    let outu_path = args.outu.unwrap_or(String::from("/dev/null"));

    let outm2_path = args.outm2.unwrap_or(String::from("/dev/null"));
    let outu2_path = args.outu2.unwrap_or(String::from("/dev/null"));

    println!(
        "Using {} threads to process reads from {}",
        num_threads, in_path
    );

    match process_reads(
        in_path,
        in2_path,
        kmer_processor,
        &outm_path,
        &outu_path,
        &outm2_path,
        &outu2_path,
        process_mode,
        ordered_output,
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

            println!(
                "\nInput:\t\t\t{} reads\t\t\t{} bases",
                read_count,
                mbase_count + ubase_count
            );
            println!(
                "Matches:\t\t{} reads ({:.2}%) \t\t{} bases ({:.2}%)",
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

    Ok(())
}

fn load_serialized_kmers(
    bin_kmers_path: &str,
    processor: &mut KmerProcessor,
) -> Result<(), Box<dyn Error>> {
    let bin_kmers_file = File::open(bin_kmers_path)?;
    let mut reader = BufReader::new(bin_kmers_file);

    processor.ref_kmers = decode_from_std_read(&mut reader, config::standard())?;

    let size_metadata = u64::MAX ^ processor.k as u64;
    if !processor.ref_kmers.contains(&size_metadata) {
        processor.ref_kmers.clear();
        return Err(format!("k-mers are of different length").into());
    }

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
    reads2_path: &Option<String>,
    matched2_path: &Option<String>,
    unmatched2_path: &Option<String>,
    interleaved_input: bool,
) -> ProcessMode {
    if reads2_path.is_some() {
        assert!(
            !interleaved_input,
            "Please disable the --interinput flag if providing 2 input files"
        );
        if matched2_path.is_none() && unmatched2_path.is_none() {
            println!(
                "Forcing interleaved output because paired input was specified for single output files"
            );
            ProcessMode::PairedInInterOut
        } else {
            assert!(
                matched2_path.is_some(),
                "Please add a second matched output path using: --outm2 <file>"
            );
            assert!(
                unmatched2_path.is_some(),
                "Please add a second unmatched output path using: --outu2 <file>"
            );
            println!("Input and output is processed as paired");
            ProcessMode::Paired
        }
    } else if interleaved_input {
        if matched2_path.is_none() && unmatched2_path.is_none() {
            println!("Input and output is processed as interleaved");
            ProcessMode::Interleaved
        } else {
            println!("Processing interleaved input and paired output");
            ProcessMode::InterInPairedOut
        }
    } else if matched2_path.is_some() || unmatched2_path.is_some() {
        panic!("Please enable the --interinput flag for 1 input file with paired output files");
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
    process_mode: ProcessMode,
    ordered_output: bool,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let processor = Arc::new(processor);

    let (chunk_sender, chunk_receiver): (
        Sender<(u32, Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)>)>,
        Receiver<(u32, Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)>)>,
    ) = channel();

    let chunk_idx = Arc::new(AtomicU32::new(0));

    let matched_filetype = matched_path.rsplit('.').next().unwrap_or("");
    let unmatched_filetype = unmatched_path.rsplit('.').next().unwrap_or("");
    let matched2_filetype = matched2_path.rsplit('.').next().unwrap_or("");
    let unmatched2_filetype = unmatched2_path.rsplit('.').next().unwrap_or("");
    let matched_stdout = matched_path == "stdout" || matched_path.starts_with("stdout.");
    let unmatched_stdout = unmatched_path == "stdout" || unmatched_path.starts_with("stdout.");
    let matched2_stdout = matched2_path == "stdout" || matched2_path.starts_with("stdout.");
    let unmatched2_stdout = unmatched2_path == "stdout" || unmatched2_path.starts_with("stdout.");

    let parallel_sender = chunk_sender.clone();
    let parellel_chunk_idx = chunk_idx.clone();
    let worker_thread = thread::spawn(move || -> Result<(), Box<dyn Error + Send + Sync>> {
        let mut reader = if reads_path == "stdin" || reads_path.starts_with("stdin.") {
            needletail::parse_fastx_reader(io::stdin())
        } else {
            needletail::parse_fastx_file(&reads_path)
        }?; // Note the trailing ? for error handling

        const CHUNK_CAP: usize = 10_000;
        let mut chunk = Vec::with_capacity(CHUNK_CAP);

        if process_mode == ProcessMode::Unpaired {
            while let Some(record) = reader.next() {
                let record = record?;
                chunk.push((
                    record.id().to_vec(),
                    record.seq().to_vec(),
                    record.qual().unwrap().to_vec(),
                ));

                if chunk.len() == CHUNK_CAP {
                    // pair idx with each chunk

                    let processor = processor.clone();
                    let sender = parallel_sender.clone();
                    let local_chunk = take(&mut chunk);

                    let current_chunk_idx = parellel_chunk_idx.fetch_add(1, Ordering::SeqCst);

                    rayon::spawn(move || {
                        let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = local_chunk
                            .into_par_iter()
                            .map(|(id, seq, qual)| {
                                let has_match = processor.process_read(&seq);
                                (has_match, id, seq, qual)
                            })
                            .collect();
                        let _ = sender.send((current_chunk_idx, processed));
                    });
                }
            }

            if !chunk.is_empty() {
                let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = chunk
                    .into_par_iter()
                    .map(|(id, seq, qual)| {
                        let has_match = processor.process_read(&seq);
                        (has_match, id, seq, qual)
                    })
                    .collect();

                let _ = parallel_sender.send((u32::MAX, processed));
            }
        } else if process_mode == ProcessMode::Interleaved
            || process_mode == ProcessMode::InterInPairedOut
        {
            while let Some(record) = reader.next() {
                let record = record?;
                chunk.push((
                    record.id().to_vec(),
                    record.seq().to_vec(),
                    record.qual().unwrap().to_vec(),
                ));

                if chunk.len() == CHUNK_CAP {
                    let processor = processor.clone();
                    let sender = parallel_sender.clone();
                    let local_chunk = take(&mut chunk);
                    let current_chunk_idx = parellel_chunk_idx.fetch_add(1, Ordering::SeqCst);

                    // Process chunk in parallel
                    rayon::spawn(move || {
                        let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = local_chunk
                            .into_par_iter()
                            .map(|(id, seq, qual)| {
                                let has_match = processor.process_read(&seq);
                                (has_match, id, seq, qual)
                            })
                            .collect();

                        let _ = sender.send((current_chunk_idx, processed));
                    });
                }
            }

            if !chunk.is_empty() {
                let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = chunk
                    .into_par_iter()
                    .map(|(id, seq, qual)| {
                        let has_match = processor.process_read(&seq);
                        (has_match, id, seq, qual)
                    })
                    .collect();

                let _ = parallel_sender.send((u32::MAX, processed));
            }
        } else if process_mode == ProcessMode::PairedInInterOut
            || process_mode == ProcessMode::Paired
        {
            let mut reader2 = parse_fastx_file(reads2_path)?;
            let mut chunk = Vec::with_capacity(CHUNK_CAP);
            while let (Some(record1), Some(record2)) = (reader.next(), reader2.next()) {
                let record1 = record1?;
                let record2 = record2?;

                chunk.push((
                    record1.id().to_vec(),
                    record1.seq().to_vec(),
                    record1.qual().unwrap().to_vec(),
                ));
                chunk.push((
                    record2.id().to_vec(),
                    record2.seq().to_vec(),
                    record2.qual().unwrap().to_vec(),
                ));

                if chunk.len() == CHUNK_CAP {
                    let processor = processor.clone();
                    let sender = parallel_sender.clone();
                    let local_chunk = take(&mut chunk);
                    let current_chunk_idx = parellel_chunk_idx.fetch_add(1, Ordering::SeqCst);

                    rayon::spawn(move || {
                        let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = local_chunk
                            .into_par_iter()
                            .map(|(id, seq, qual)| {
                                let has_match = processor.process_read(&seq);
                                (has_match, id, seq, qual)
                            })
                            .collect();

                        let _ = sender.send((current_chunk_idx, processed));
                    });
                }
            }

            if !chunk.is_empty() {
                let processed: Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)> = chunk
                    .into_par_iter()
                    .map(|(id, seq, qual)| {
                        let has_match = processor.process_read(&seq);
                        (has_match, id, seq, qual)
                    })
                    .collect();

                let _ = parallel_sender.send((u32::MAX, processed));
            }
        }

        Ok(())
    });

    drop(chunk_sender);

    let mut matched_writer: BufWriter<File> = BufWriter::new(File::create(matched_path)?);
    let mut unmatched_writer: BufWriter<File> =
        BufWriter::with_capacity(4_000_000, File::create(unmatched_path)?);

    let matched_count = Arc::new(AtomicU32::new(0));
    let matched_bases = Arc::new(AtomicU32::new(0));

    let unmatched_count = Arc::new(AtomicU32::new(0));
    let unmatched_bases = Arc::new(AtomicU32::new(0));

    let m_count;
    let m_bases;
    let u_count;
    let u_bases;

    if ordered_output {
        let mut output_chunks = Vec::new();
        for output_chunk in chunk_receiver {
            output_chunks.push(output_chunk);
        }

        output_chunks.sort_by(|a, b| a.0.cmp(&b.0));

        (m_count, m_bases, u_count, u_bases) = process_output_chunks(
            output_chunks.into_iter(),
            &mut matched_writer,
            matched_count,
            matched_bases,
            matched_filetype,
            matched_stdout,
            matched2_path,
            matched2_filetype,
            matched2_stdout,
            &mut unmatched_writer,
            unmatched_count,
            unmatched_bases,
            unmatched_filetype,
            unmatched_stdout,
            unmatched2_path,
            unmatched2_filetype,
            unmatched2_stdout,
            process_mode,
        )?;
    } else {
        (m_count, m_bases, u_count, u_bases) = process_output_chunks(
            chunk_receiver.into_iter(),
            &mut matched_writer,
            matched_count,
            matched_bases,
            matched_filetype,
            matched_stdout,
            matched2_path,
            matched2_filetype,
            matched2_stdout,
            &mut unmatched_writer,
            unmatched_count,
            unmatched_bases,
            unmatched_filetype,
            unmatched_stdout,
            unmatched2_path,
            unmatched2_filetype,
            unmatched_stdout,
            process_mode,
        )?;
    }

    matched_writer.flush()?;
    unmatched_writer.flush()?;
    worker_thread.join().unwrap()?;

    // delete file is stdout is input
    for path in [
        &matched_path,
        unmatched_path,
        matched2_path,
        unmatched2_path,
    ] {
        if path.starts_with("stdout") {
            if let Err(e) = fs::remove_file(path) {
                // Ignore error if file doesn't exist
                if e.kind() != io::ErrorKind::NotFound {
                    eprintln!("Warning: Could not delete file '{}': {}", path, e);
                }
            }
        }
    }

    Ok((m_count, m_bases, u_count, u_bases))
}

fn process_output_chunks(
    chunks_iterator: impl Iterator<Item = (u32, Vec<(bool, Vec<u8>, Vec<u8>, Vec<u8>)>)>,
    m_writer: &mut BufWriter<File>,
    m_count: Arc<AtomicU32>,
    m_bases: Arc<AtomicU32>,
    m_filetype: &str,
    m_stdout: bool,
    m2_path: &str,
    m2_filetype: &str,
    m2_stdout: bool,
    u_writer: &mut BufWriter<File>,
    u_count: Arc<AtomicU32>,
    u_bases: Arc<AtomicU32>,
    u_filetype: &str,
    u_stdout: bool,
    u2_path: &str,
    u2_filetype: &str,
    u2_stdout: bool,
    process_mode: ProcessMode,
) -> Result<(u32, u32, u32, u32), Box<dyn Send + Sync + Error>> {
    if process_mode == ProcessMode::Unpaired {
        for (_, chunk) in chunks_iterator {
            for (has_match, id, seq, qual) in chunk {
                if has_match {
                    write_read(m_writer, &id, &seq, &qual, m_filetype, m_stdout)?;

                    m_count.fetch_add(1, Ordering::Relaxed);
                    m_bases.fetch_add(seq.len() as u32, Ordering::Relaxed);
                } else {
                    write_read(u_writer, &id, &seq, &qual, u_filetype, u_stdout)?;

                    u_count.fetch_add(1, Ordering::Relaxed);
                    u_bases.fetch_add(seq.len() as u32, Ordering::Relaxed);
                }
            }
        }
    } else if process_mode == ProcessMode::Interleaved
        || process_mode == ProcessMode::PairedInInterOut
    {
        for (_, chunk) in chunks_iterator {
            for i in (0..chunk.len() - 1).step_by(2) {
                let has_match = chunk[i].0 || chunk[i + 1].0;

                if has_match {
                    write_read(
                        m_writer,
                        &chunk[i].1,
                        &chunk[i].2,
                        &chunk[i].3,
                        m_filetype,
                        m_stdout,
                    )?;

                    write_read(
                        m_writer,
                        &chunk[i + 1].1,
                        &chunk[i + 1].2,
                        &chunk[i + 1].3,
                        m_filetype,
                        m_stdout,
                    )?;

                    m_count.fetch_add(2, Ordering::Relaxed);
                    m_bases.fetch_add(
                        (chunk[i].2.len() + chunk[i + 1].2.len()) as u32,
                        Ordering::Relaxed,
                    );
                } else {
                    write_read(
                        u_writer,
                        &chunk[i].1,
                        &chunk[i].2,
                        &chunk[i].3,
                        u_filetype,
                        u_stdout,
                    )?;

                    write_read(
                        u_writer,
                        &chunk[i + 1].1,
                        &chunk[i + 1].2,
                        &chunk[i + 1].3,
                        u_filetype,
                        u_stdout,
                    )?;

                    u_count.fetch_add(2, Ordering::Relaxed);
                    u_bases.fetch_add(
                        (chunk[i].2.len() + chunk[i + 1].2.len()) as u32,
                        Ordering::Relaxed,
                    );
                }
            }
        }
    } else if process_mode == ProcessMode::InterInPairedOut || process_mode == ProcessMode::Paired {
        let mut m2_writer: BufWriter<File> = BufWriter::new(File::create(m2_path)?);
        let mut u2_writer: BufWriter<File> = BufWriter::new(File::create(u2_path)?);

        for (_, chunk) in chunks_iterator {
            for i in (0..chunk.len() - 1).step_by(2) {
                let has_match = chunk[i].0 || chunk[i + 1].0;

                if has_match {
                    write_read(
                        m_writer,
                        &chunk[i].1,
                        &chunk[i].2,
                        &chunk[i].3,
                        m_filetype,
                        m_stdout,
                    )?;

                    write_read(
                        &mut m2_writer,
                        &chunk[i + 1].1,
                        &chunk[i + 1].2,
                        &chunk[i + 1].3,
                        m2_filetype,
                        m2_stdout,
                    )?;

                    m_count.fetch_add(2, Ordering::Relaxed);
                    m_bases.fetch_add(
                        (chunk[i].2.len() + chunk[i + 1].2.len()) as u32,
                        Ordering::Relaxed,
                    );
                } else {
                    write_read(
                        u_writer,
                        &chunk[i].1,
                        &chunk[i].2,
                        &chunk[i].3,
                        u_filetype,
                        u_stdout,
                    )?;

                    write_read(
                        &mut u2_writer,
                        &chunk[i + 1].1,
                        &chunk[i + 1].2,
                        &chunk[i + 1].3,
                        u2_filetype,
                        u2_stdout,
                    )?;

                    u_count.fetch_add(2, Ordering::Relaxed);
                    u_bases.fetch_add(
                        (chunk[i].2.len() + chunk[i + 1].2.len()) as u32,
                        Ordering::Relaxed,
                    );
                }
            }
        }

        m2_writer.flush()?;
        u2_writer.flush()?;
    }

    Ok((
        m_count.load(Ordering::Relaxed),
        m_bases.load(Ordering::Relaxed),
        u_count.load(Ordering::Relaxed),
        u_bases.load(Ordering::Relaxed),
    ))
}

fn write_read(
    writer: &mut BufWriter<File>,
    id: &[u8],
    sequence: &[u8],
    quality: &[u8],
    format: &str,
    stdout: bool,
) -> Result<(), Box<dyn Send + Sync + Error>> {
    if stdout {
        unsafe {
            if format == "fa" || format == "fna" || format == "fasta" {
                println!(">");
                println!("{}", str::from_utf8_unchecked(id));
                println!("{}", str::from_utf8_unchecked(sequence));
            } else {
                println!(">");
                println!("{}", str::from_utf8_unchecked(id));
                println!("{}\n+", str::from_utf8_unchecked(sequence));
                println!("{}", str::from_utf8_unchecked(quality));
            }
        }
    } else {
        if format == "fa" || format == "fna" || format == "fasta" {
            writer.write_all(b">")?;
            writer.write_all(id)?;
            writer.write_all(b"\n")?;
            writer.write_all(sequence)?;
            writer.write_all(b"\n")?;
        } else {
            writer.write_all(b"@")?;
            writer.write_all(id)?;
            writer.write_all(b"\n")?;
            writer.write_all(sequence)?;
            writer.write_all(b"\n+\n")?;
            writer.write_all(quality)?;
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}
