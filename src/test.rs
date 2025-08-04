use crate::kmer_processor::*;
use rand::Rng;
use rustc_hash::FxHashSet;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

pub fn reverse_complement_str(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => 'X',
        })
        .collect()
}

#[test]
fn make_test_file() -> Result<(), Box<dyn std::error::Error>> {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut rand = rand::rng();
    let mut seq = String::new();

    let matched_file = File::create("in/1m_150.fq")?;
    let mut matched_writer = BufWriter::new(matched_file);

    for x in 0..1_000_000 {
        for _i in 0..150 {
            let base = BASES[rand.random_range(0..3)];
            seq.push(base);
        }

        writeln!(matched_writer, "@{}", x)?;
        writeln!(matched_writer, "{}", seq)?;
        writeln!(matched_writer, "+")?;
        writeln!(
            matched_writer,
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        )?;
        seq = String::new();
    }

    Ok(())
}

// TESTS
// #[test]
// fn test_kmer_fns() {
//     let seq_vec = b"TGCTCAGATCATGTTTGTGTGG";
//     let kmer = encode(seq_vec[0..21].try_into().unwrap());
//     let kmer2 = encode(seq_vec[1..22].try_into().unwrap());

//     println!("Encoded k-mer: {:042b}", kmer);
//     println!("Encoded k-mer: {:042b}", kmer2);
//     assert_ne!(kmer, kmer2);
//     assert_eq!(kmer, 0b111001110100100011010011101111111011101110);

//     let bit_cap = (1u64 << 21 * 2) - 1;
//     let mut shifted_kmer = kmer;
//     assert_eq!(shifted_kmer, kmer);

//     shifted_kmer = ((kmer << 2) | encode(b"G")) & bit_cap;
//     println!("Shifted k-mer: {:b}", shifted_kmer);

//     assert_eq!(shifted_kmer, kmer2);
//     assert_ne!(shifted_kmer, kmer);

//     println!("Sequence vector: {:?}", seq_vec);

//     let decoded_seq = [
//         84, 71, 67, 84, 67, 65, 71, 65, 84, 67, 65, 84, 71, 84, 84, 84, 71, 84, 71, 84, 71, 71,
//     ];

//     assert_eq!(decode(kmer, 21), decoded_seq[0..21]);
//     println!("Decoded k-mer 1: {:?}", decode(kmer, 21));

//     assert_eq!(decode(kmer2, 21), decoded_seq[1..22]);
//     println!("Decoded k-mer 2: {:?}", decode(kmer2, 21));

//     assert_eq!(decode(shifted_kmer, 21), decoded_seq[1..22]);
//     println!("Decoded shifter: {:?}", decode(shifted_kmer, 21));
// }

// #[test]
// fn test_kmer_processor() {
//     let k = 21;
//     let threshold = 1;
//     let mut processor = KmerProcessor::new(k, threshold);
//     let mut test_hash = FxHashSet::default();

//     // Test reference processing
//     let ref_seq = b"TGCTCAGATCATGTTTGTGTGAA";
//     let rev_seq = b"TTCACACAAACATGATCTGAGCA";

//     processor.process_ref(&ref_seq.to_vec());

//     test_hash.insert(canonical_kmer(encode(&ref_seq[0..k]), k));
//     println!("First encoded K-mer: {}", encode(&ref_seq[0..k]));
//     test_hash.insert(canonical_kmer(encode(&ref_seq[1..k + 1]), k));
//     println!("Second encoded K-mer: {}", encode(&ref_seq[1..k + 1]));
//     test_hash.insert(canonical_kmer(encode(&ref_seq[2..k + 2]), k));

//     assert!(processor.process_read(&ref_seq.to_vec()));
//     assert!(processor.process_read(&rev_seq.to_vec()));

//     println!("Ref K-mers in KP: {:#?}", processor.ref_kmers);
//     println!("Test hash K-mers: {:#?}", test_hash);
//     assert_eq!(processor.ref_kmers, test_hash);

//     println!("Original k-mer: {:b}", encode(&ref_seq[0..k]));
//     println!(
//         "Canon OG k-mer: {:b}",
//         canonical_kmer(encode(&ref_seq[0..k]), k)
//     );
//     println!("RC k-mer:       {:b}", encode(&rev_seq[2..k + 2]));
//     println!(
//         "Canon RC k-mer: {:b}",
//         canonical_kmer(encode(&rev_seq[2..k + 2]), k)
//     );

//     assert_eq!(
//         encode(&rev_seq[2..k + 2]),
//         canonical_kmer(encode(&rev_seq[2..k + 2]), k)
//     );
//     assert_eq!(
//         encode(&rev_seq[2..k + 2]),
//         canonical_kmer(encode(&ref_seq[0..k]), k)
//     );
//     assert_eq!(
//         canonical_kmer(encode(&rev_seq[0..k]), k),
//         canonical_kmer(encode(&ref_seq[2..k + 2]), k)
//     );

//     assert!(processor.ref_kmers.contains(&canonical_kmer(
//         encode(ref_seq[0..k].try_into().unwrap()),
//         k
//     )));

//     // Test read processing
//     let read_seq = b"TGCTCAGATCATGTTTGTGTGG";
//     assert!(processor.process_read(read_seq));

//     // Test read with insufficient k-mers
//     let short_read = b"TGCTCAGATC";
//     assert!(!processor.process_read(short_read));
// }

// #[test]
// fn test_misc() {
//     let seq = "CCGACG";
//     let rev = reverse_complement_str(&reverse_complement_str(seq));
//     println!("Original sequence:  {}\nReverse complement: {}", seq, rev);
//     assert_eq!(seq, rev);

//     let arr_seq = b"AAAAA";
//     let processed_seq = decode(encode(arr_seq), 5);
//     println!("Original: {:#?}\nProcessed: {:#?}", arr_seq, processed_seq);
//     assert_eq!(arr_seq.to_vec(), processed_seq);
// }

// pub fn decode(encoded: u64, k: usize) -> Vec<u8> {
//     let mut seq = Vec::new();
//     for i in 0..k {
//         let base = match (encoded >> (i * 2)) & 0b11 {
//             0b00 => b'A',
//             0b01 => b'C',
//             0b10 => b'G',
//             0b11 => b'T',
//             _ => panic!("Non-DNA base!"),
//         };
//         seq.push(base);
//     }
//     seq.into_iter().rev().collect()
// }

/* OLD CODE
match load_kmer_index(serialized_kmers_filename, &mut kmer_processor) {
    Ok(kmers) => {
        println!(
            "Loaded {} k-mers from pre-built index file: {}",
            kmer_processor.ref_kmers.len(),
            serialized_kmers_filename
        );
        kmers
    }
    Err(_) => {
        println!("No pre-built index found, building from reference sequences...");

        // Load reference sequences (can use streaming for large reference files)
        let ref_seqs = match load_reference_sequences(reference_filename, &mut kmer_processor) {
            Ok(_ok) => {
                println!("Loaded reference sequences");
            }
            Err(e) => {
                eprintln!("Error loading reference sequences: {}", e);
                std::process::exit(1);
            }
        };

        // Save for future use
        if let Err(e) =
            save_kmer_index(kmer_processor.ref_kmers.clone(), serialized_kmers_filename)
        {
            eprintln!("Warning: Could not save k-mer index: {}", e);
        } else {
            println!("Saved k-mer index for future use");
        }

        ref_seqs
    }
};

Process reads using streaming (much more memory efficient)
match process_reads(&query_filename, kmer_processor) {
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

fn process_reads(
    reads_path: &str,
    processor: KmerProcessor,
) -> Result<(Vec<(String, String)>, Vec<(String, String)>), Box<dyn std::error::Error>> {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();

    // needletail handles file format detection and streaming
    let mut reader = needletail::parse_fastx_file(reads_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        let seq = String::from_utf8_lossy(&record.seq()).to_string();

        let is_match = processor.process_read(&seq);

        // Classify this read based on threshold
        if is_match {
            matched.push((id, seq));
        } else {
            unmatched.push((id, seq));
        }
    }

    Ok((matched, unmatched))
}

#[inline(always)]
pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut shift = (k - 1) * 2;
    for i in 0..k {
        let base = (kmer >> (i * 2)) & 0b11;
        let comp = base ^ 0b11;
        rc |= comp << shift;
        shift -= 2;
    }
    rc
}

#[inline(always)]
pub fn canonical_kmer(kmer: u64, k: usize) -> u64 {
    let rc = reverse_complement(kmer, k);
    if kmer < rc { kmer } else { rc }
}

pub fn encode_old(sequence: &[u8]) -> u64 {
    sequence.iter().fold(0, |acc, &base| {
        (acc << 2)
            | match base {
                b'A' => 0b00,
                b'C' => 0b01,
                b'G' => 0b10,
                b'T' => 0b11,
                _ => panic!("Invalid nucleotide base"),
            }
    })
}

fn process_reads_thread(
    reads_path: &str,
    processor: KmerProcessor,
    num_threads: usize,
    matched_path: &str,
    unmatched_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let processor = Arc::from(processor);

    let matched: Arc<Vec<&str>> = Arc::new(Vec::new());
    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_path)?)));
    let mut matched_count = Arc::new(Mutex::new(0u32));

    let unmatched: Arc<Vec<&str>> = Arc::new(Vec::new());
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_path)?)));
    let mut unmatched_count = Arc::new(Mutex::new(0u32));

    // Create a shared queue for distributing reads to worker threads
    let queue: Arc<Mutex<Vec<(Arc<str>, Arc<[u8]>)>>> = Arc::new(Mutex::new(Vec::new()));

    // Read and add sequences to queue using needletail
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    while let Some(record) = reader.next() {
        let record = record?;
        let id: Arc<str> = Arc::from(str::from_utf8(record.id())?);
        let seq: Arc<[u8]> = Arc::from(record.seq());
        queue.lock().unwrap().push((id, seq));
    }

    // Spawn worker threads
    let mut handles = Vec::new();
    for _ in 0..num_threads {
        let queue = Arc::clone(&queue);
        let processor = Arc::clone(&processor);
        let matched_writer = Arc::clone(&matched_writer);
        let unmatched_writer = Arc::clone(&unmatched_writer);
        let matched_count = Arc::clone(&matched_count);
        let unmatched_count = Arc::clone(&unmatched_count);

        let handle = thread::spawn(move || {
            loop {
                let read = {
                    let mut queue = queue.lock().unwrap();
                    queue.pop()
                };

                if let Some(read) = read {
                    let id = read.0;
                    let seq = read.1;
                    let is_match = processor.process_read(&seq);
                    let mut matched_count = matched_count.lock().unwrap();
                    let mut unmatched_count = unmatched_count.lock().unwrap();

                    if is_match {
                        let mut writer = matched_writer.lock().unwrap();
                        writeln!(writer, ">{}", id).unwrap();
                        writer.write_all(&seq).unwrap();
                        writeln!(writer).unwrap();
                        *matched_count += 1;
                    } else {
                        let mut writer = unmatched_writer.lock().unwrap();
                        writeln!(writer, ">{}", id).unwrap();
                        writer.write_all(&seq).unwrap();
                        writeln!(writer).unwrap();
                        *unmatched_count += 1;
                    }
                } else {
                    break;
                }
            }
        });

        handles.push(handle);
    }

    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }

    println!("{}", matched_count.lock().unwrap());
    println!("{}", unmatched_count.lock().unwrap());

    matched_writer.lock().unwrap().flush()?;
    unmatched_writer.lock().unwrap().flush()?;

    Ok(())
}

fn process_reads(
    reads_path: &str,
    processor: KmerProcessor,
    matched_path: &str,
    unmatched_path: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let chunk_size = 100_000;
    let chunks = get_chunks(reads_path, chunk_size)?;

    let processor = Arc::new(processor);
    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_path)?)));
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_path)?)));

    let mseq_count = Arc::new(Mutex::new(0u32));
    let mbase_count = Arc::new(Mutex::new(0u32));

    let useq_count = Arc::new(Mutex::new(0u32));
    let ubase_count = Arc::new(Mutex::new(0u32));

    chunks.into_par_iter().try_for_each(|chunk| {
        for record in &chunk.reads {
            let seq = &chunk.buffer[record.range.clone()];
            let is_match = processor.process_read(seq);

            if is_match {
                let mut writer = matched_writer.lock().unwrap();
                writeln!(writer, ">{}", record.id)?;
                writer.write_all(seq)?;
                writeln!(writer)?;
                *mseq_count.lock().unwrap() += 1;
                *mbase_count.lock().unwrap() += seq.len() as u32;
            } else {
                let mut writer = unmatched_writer.lock().unwrap();
                writeln!(writer, ">{}", record.id)?;
                writer.write_all(seq)?;
                writeln!(writer)?;
                *useq_count.lock().unwrap() += 1;
                *ubase_count.lock().unwrap() += seq.len() as u32;
            }
        }
        Ok::<(), Box<dyn Error + Send + Sync>>(())
    })?;

    let mseq_count = *mseq_count.lock().unwrap();
    let mbase_count = *mbase_count.lock().unwrap();

    let useq_count = *useq_count.lock().unwrap();
    let ubase_count = *ubase_count.lock().unwrap();

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn get_chunks(
    reads_path: &str,
    chunk_size: usize,
) -> Result<Vec<Chunk>, Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;

    let mut chunks = Vec::new();
    let mut reads = Vec::with_capacity(chunk_size);
    let mut buffer = Vec::with_capacity(chunk_size * 150);

    while let Some(record) = reader.next() {
        let record = record?;
        let id = Arc::from(str::from_utf8(record.id())?);
        let seq = record.seq();

        let start = buffer.len();
        buffer.extend_from_slice(&seq);
        let end = buffer.len();

        reads.push(RecordRef::new(id, start..end));

        if reads.len() == chunk_size {
            chunks.push(Chunk::new(buffer, reads));
            reads = Vec::with_capacity(chunk_size);
            buffer = Vec::with_capacity(chunk_size * 150);
        }
    }

    if !reads.is_empty() {
        chunks.push(Chunk::new(buffer, reads));
    }

    Ok(chunks)
}

*/
