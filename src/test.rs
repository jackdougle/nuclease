use bio::pattern_matching::shift_and::ShiftAnd;
use rustc_hash::FxHashSet;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};

use crate::kmer::*;
use crate::kmer_processor::KmerProcessor;

pub fn run() {
    let chars: Vec<char> = vec!['A', 'T', 'C', 'G'];
    let rng = 2;

    for _x in 0..10 {
        let mut word = String::new();
        // Use rand::Rng's gen_range method and a proper rng instance
        while word.len() <= 120 {
            let n = rng;
            word.push(chars[n]);
        }
        println!("{}", word);
    }
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

fn extract_kmers(seq: &str, k: usize, canonical: bool) -> HashSet<String> {
    let mut kmers = HashSet::new();
    for i in 0..=seq.len() - k {
        let kmer = &seq[i..i + k];
        if canonical {
            let rc = reverse_complement(kmer);
            kmers.insert(std::cmp::min(kmer, &rc).to_string());
        } else {
            kmers.insert(kmer.to_string());
        }
    }
    kmers
}

fn load_kmers(path: &str, k: usize, canonical: bool) -> HashSet<String> {
    let file = File::open(path).expect("Cannot open reference file");
    let reader = BufReader::new(file);

    let mut kmers = HashSet::new();
    let mut seq = String::new();
    for line in reader.lines() {
        let l = line.unwrap();
        if l.starts_with('>') {
            if !seq.is_empty() {
                kmers.extend(extract_kmers(&seq, k, canonical));
                seq.clear();
            }
        } else {
            seq.push_str(&l);
        }
    }
    if !seq.is_empty() {
        kmers.extend(extract_kmers(&seq, k, canonical));
    }

    kmers
}

fn classify_reads(
    reads: Vec<String>,
    k: usize,
    kmers: &HashSet<String>,
    threshold: usize,
    canonical: bool,
) -> (Vec<String>, Vec<String>) {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();

    for read in reads {
        let read_kmers = extract_kmers(&read, k, canonical);
        let hits = read_kmers.intersection(kmers).count();
        if hits >= threshold {
            matched.push(read);
        } else {
            unmatched.push(read);
        }
    }

    (matched, unmatched)
}

pub fn main() {
    let reference_path = "input/refs.fasta";
    let reads = vec![
        "GACGATCGATCGATCGATTACTTACTCAGGACGCTAG".to_string(),
        "AGCTAGCTAGCTACGTACGCACTGCTAGCTAGCTAGCT".to_string(),
    ];

    // let test_path = "samples.fastq";

    let k = 21;
    let threshold = 3;
    let canonical = true;

    let reference_kmers = load_kmers(reference_path, k, true);
    println!("Reference kmers: {:?}", reference_kmers);

    // let test_kmers: HashSet<String> = load_kmers(test_path, k, false);
    // println!("Test kmers: {:?}", test_kmers);

    let (matched, unmatched) = classify_reads(reads, k, &reference_kmers, threshold, canonical);

    println!("Matched reads: {:?}", matched);
    println!("Unmatched reads: {:?}", unmatched);
}

fn get_read_kmers(read_seqs: &HashMap<String, String>, k: usize) -> HashMap<String, Vec<String>> {
    let mut read_kmers: HashMap<String, Vec<String>> = HashMap::new();

    for (name, seq) in read_seqs {
        let expected_kmers = seq.len().saturating_sub(k) + 1;
        let kmers_vec = read_kmers
            .entry(name.clone())
            .or_insert_with(|| Vec::with_capacity(expected_kmers));

        for i in 0..=seq.len() - k {
            let kmer = &seq[i..i + k];

            kmers_vec.push(kmer.to_string());
        }
    }

    read_kmers
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
                            let rc = reverse_complement(kmer);
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
    println!("First encoded K-mer: {}", encode(&ref_seq[1..k + 1]));
    test_hash.insert(canonical_kmer(encode(&ref_seq[2..k + 2]), k));

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
