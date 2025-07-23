use crate::kmer_processor::*;
use bio::pattern_matching::shift_and::ShiftAnd;
use rustc_hash::FxHashSet;
use std::collections::HashSet;
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

/// Processes reads with parallel threading for better performance
///
/// PURPOSE: Implement your goal of "number of CPUs to use" with "at least some
/// parallelism (say up to 4cpu) is essential". This function demonstrates how
/// to use needletail with threading for parallel processing.
///
/// HOW IT WORKS WITH NEEDLETAIL:
/// - Uses needletail's streaming to read reads efficiently
/// - Distributes reads across multiple worker threads
/// - Each thread processes its assigned reads independently
/// - Combines results from all threads
///
/// BENEFITS:
/// - Utilizes multiple CPU cores
/// - Maintains memory efficiency through streaming
/// - Scales with available hardware
/// - Supports your goal of essential parallelism
fn process_reads_parallel(
    reads_path: &str,
    ref_kmers: Arc<HashSet<String>>,
    k: usize,
    threshold: usize,
    canonical: bool,
    num_threads: usize,
) -> Result<(Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    use std::sync::{Arc, Mutex};
    use std::thread;

    let matched = Arc::new(Mutex::new(Vec::new()));
    let unmatched = Arc::new(Mutex::new(Vec::new()));

    // Create a shared queue for distributing reads to worker threads
    let queue: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let queue_clone = Arc::clone(&queue);

    // Spawn worker threads
    let mut handles = Vec::new();
    for _ in 0..num_threads {
        let queue = Arc::clone(&queue);
        let ref_kmers = Arc::clone(&ref_kmers);
        let matched = Arc::clone(&matched);
        let unmatched = Arc::clone(&unmatched);

        let handle = thread::spawn(move || {
            loop {
                let seq = {
                    let mut queue = queue.lock().unwrap();
                    queue.pop()
                };

                if let Some(seq) = seq {
                    let mut read_kmers = HashSet::new();
                    for i in 0..=seq.len() - k {
                        let kmer = &seq[i..i + k];

                        if canonical {
                            let rc = reverse_complement_str(kmer);
                            read_kmers.insert(std::cmp::min(kmer, &rc).to_string());
                        } else {
                            read_kmers.insert(kmer.to_string());
                        }
                    }

                    let hits = read_kmers.intersection(&ref_kmers).count();

                    if hits >= threshold {
                        matched.lock().unwrap().push(seq);
                    } else {
                        unmatched.lock().unwrap().push(seq);
                    }
                } else {
                    break; // No more work to do
                }
            }
        });

        handles.push(handle);
    }

    // Read and add sequences to queue using needletail
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = String::from_utf8_lossy(&record.seq()).to_string();
        queue_clone.lock().unwrap().push(seq);
    }

    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }

    let matched = Arc::try_unwrap(matched).unwrap().into_inner().unwrap();
    let unmatched = Arc::try_unwrap(unmatched).unwrap().into_inner().unwrap();

    Ok((matched, unmatched))
}

#[test]
fn test_kmer_fns() {
    let seq_vec = b"TGCTCAGATCATGTTTGTGTGG";
    let kmer = encode(seq_vec[0..21].try_into().unwrap());
    let kmer2 = encode(seq_vec[1..22].try_into().unwrap());

    println!("Encoded k-mer: {:042b}", kmer);
    println!("Encoded k-mer: {:042b}", kmer2);
    assert_ne!(kmer, kmer2);
    assert_eq!(kmer, 0b111001110100100011010011101111111011101110);

    let bit_cap = (1u64 << 21 * 2) - 1;
    let mut shifted_kmer = kmer;
    assert_eq!(shifted_kmer, kmer);

    shifted_kmer = ((kmer << 2) | encode(b"G")) & bit_cap;
    println!("Shifted k-mer: {:b}", shifted_kmer);

    assert_eq!(shifted_kmer, kmer2);
    assert_ne!(shifted_kmer, kmer);

    println!("Sequence vector: {:?}", seq_vec);

    let decoded_seq = [
        84, 71, 67, 84, 67, 65, 71, 65, 84, 67, 65, 84, 71, 84, 84, 84, 71, 84, 71, 84, 71, 71,
    ];

    assert_eq!(decode(kmer, 21), decoded_seq[0..21]);
    println!("Decoded k-mer 1: {:?}", decode(kmer, 21));

    assert_eq!(decode(kmer2, 21), decoded_seq[1..22]);
    println!("Decoded k-mer 2: {:?}", decode(kmer2, 21));

    assert_eq!(decode(shifted_kmer, 21), decoded_seq[1..22]);
    println!("Decoded shifter: {:?}", decode(shifted_kmer, 21));
}

#[test]
fn test_kmer_processor() {
    let k = 21;
    let threshold = 1;
    let mut processor = KmerProcessor::new(k, threshold);
    let mut test_hash = FxHashSet::default();

    // Test reference processing
    let ref_seq = b"TGCTCAGATCATGTTTGTGTGAA";
    let rev_seq = b"TTCACACAAACATGATCTGAGCA";

    processor.process_ref(ref_seq.to_vec());

    test_hash.insert(canonical_kmer(encode(&ref_seq[0..k]), k));
    println!("First encoded K-mer: {}", encode(&ref_seq[0..k]));
    test_hash.insert(canonical_kmer(encode(&ref_seq[1..k + 1]), k));
    println!("Second encoded K-mer: {}", encode(&ref_seq[1..k + 1]));
    test_hash.insert(canonical_kmer(encode(&ref_seq[2..k + 2]), k));

    assert!(processor.process_read(&String::from_utf8(ref_seq.to_vec()).unwrap()));
    assert!(processor.process_read(&String::from_utf8(rev_seq.to_vec()).unwrap()));

    println!("Ref K-mers in KP: {:#?}", processor.ref_kmers);
    println!("Test hash K-mers: {:#?}", test_hash);
    assert_eq!(processor.ref_kmers, test_hash);

    println!("Original k-mer: {:b}", encode(&ref_seq[0..k]));
    println!(
        "Canon OG k-mer: {:b}",
        canonical_kmer(encode(&ref_seq[0..k]), k)
    );
    println!("RC k-mer:       {:b}", encode(&rev_seq[2..k + 2]));
    println!(
        "Canon RC k-mer: {:b}",
        canonical_kmer(encode(&rev_seq[2..k + 2]), k)
    );

    assert_eq!(
        encode(&rev_seq[2..k + 2]),
        canonical_kmer(encode(&rev_seq[2..k + 2]), k)
    );
    assert_eq!(
        encode(&rev_seq[2..k + 2]),
        canonical_kmer(encode(&ref_seq[0..k]), k)
    );
    assert_eq!(
        canonical_kmer(encode(&rev_seq[0..k]), k),
        canonical_kmer(encode(&ref_seq[2..k + 2]), k)
    );

    assert!(processor.ref_kmers.contains(&canonical_kmer(
        encode(ref_seq[0..k].try_into().unwrap()),
        k
    )));

    // Test read processing
    let read_seq = "TGCTCAGATCATGTTTGTGTGG";
    assert!(processor.process_read(read_seq));

    // Test read with insufficient k-mers
    let short_read = "TGCTCAGATC";
    assert!(!processor.process_read(short_read));
}

#[test]
fn test_misc() {
    let seq = "CCGACG";
    let rev = reverse_complement_str(&reverse_complement_str(seq));
    println!("Original sequence:  {}\nReverse complement: {}", seq, rev);
    assert_eq!(seq, rev);

    let arr_seq = b"AAAAA";
    let processed_seq = decode(encode(arr_seq), 5);
    println!("Original: {:#?}\nProcessed: {:#?}", arr_seq, processed_seq);
    assert_eq!(arr_seq.to_vec(), processed_seq);
}

#[test]
fn test_bio() {
    let pattern = b"AGCTAGCTACG";
    let shifter = ShiftAnd::new(pattern);
}
