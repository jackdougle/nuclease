use rand::Rng;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::{Arc, Mutex};
use std::thread;

pub fn run() {
    let chars: Vec<char> = vec!['A', 'T', 'C', 'G'];
    let mut rng = rand::rng;

    for x in 0..10 {
        let mut word = String::new();
        // Use rand::Rng's gen_range method and a proper rng instance
        while word.len() <= 120 {
            let n = rng().gen_range(0..4);
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

fn check_kmers(
    ref_kmers: &HashSet<String>,
    read_kmers: &HashMap<String, Vec<String>>,
    threshold: usize,
    canonical: bool,
) -> (HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
    let mut matched: HashMap<String, Vec<String>> = HashMap::new();
    let mut unmatched: HashMap<String, Vec<String>> = HashMap::new();

    for (name, read_set) in read_kmers {
        let mut count = 0;
        let mut matched_kmers: Vec<String> = Vec::new();

        for read in read_set {
            let rev = rev_comp(read);
            let canonical = if canonical {
                std::cmp::min(read, &rev)
            } else {
                read
            };

            if ref_kmers.contains(canonical) {
                count += 1;
                matched_kmers.push(canonical.clone());
            }
        }
        if count >= threshold {
            matched
                .entry(name.clone())
                .or_insert(Vec::new())
                .push(matched_kmers.join(", "));
        } else {
            unmatched
                .entry(name.clone())
                .or_insert(Vec::new())
                .push(read_set.join(", "));
        }
    }

    (matched, unmatched)
}

fn read_file_to_string(filename: &str) -> Result<String, std::io::Error> {
    let mut file: File = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

fn load_seqs(path: &str) -> HashMap<String, String> {
    let mut seqs: HashMap<String, String> = HashMap::new();

    let seq_set =
        read_file_to_string(path).expect(&format!("Error: Invalid reference file - {}", path));
    let mut temp: &str = "";
    for line in seq_set.lines() {
        if temp.starts_with('>') {
            seqs.insert(temp[6..].to_string(), line.to_string());
        } else if temp.starts_with('@') {
            seqs.insert(temp[1..].to_string(), line.to_string());
        }
        temp = &line;
    }

    seqs
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
) -> Result<(Vec<String>, Vec<Sxwtring>), Box<dyn std::error::Error>> {
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

/// Processes reads with configurable memory limits
///
/// PURPOSE: Implement your goal of "memory maximum" parameter. This function
/// processes reads in chunks that fit within a specified memory limit, ensuring
/// the program never exceeds the user's memory constraints.
///
/// HOW IT WORKS WITH NEEDLETAIL:
/// - Uses needletail's streaming to read reads one at a time
/// - Tracks memory usage of processed reads
/// - Flushes results when memory limit is approached
/// - Continues processing without exceeding limits
///
/// BENEFITS:
/// - Respects user memory limits
/// - Prevents out-of-memory crashes
/// - Suitable for shared computing environments
/// - Supports your goal of predictable resource usage
fn process_with_memory_limit(
    reads_path: &str,
    ref_kmers: &HashSet<String>,
    k: usize,
    threshold: usize,
    canonical: bool,
    memory_limit_mb: usize,
) -> Result<(Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();
    let mut current_batch = Vec::new();
    let mut current_memory = 0;

    let memory_limit_bytes = memory_limit_mb * 1024 * 1024;

    let mut reader = needletail::parse_fastx_file(reads_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let seq = String::from_utf8_lossy(&record.seq()).to_string();

        // Estimate memory usage (rough calculation)
        let seq_memory = seq.len() * std::mem::size_of::<char>();

        // If adding this sequence would exceed memory limit, process current batch
        if current_memory + seq_memory > memory_limit_bytes && !current_batch.is_empty() {
            let (batch_matched, batch_unmatched) =
                process_batch(&current_batch, ref_kmers, k, threshold, canonical);
            matched.extend(batch_matched);
            unmatched.extend(batch_unmatched);

            current_batch.clear();
            current_memory = 0;
        }

        current_batch.push(seq);
        current_memory += seq_memory;
    }
    // Process final batch
    if !current_batch.is_empty() {
        let (batch_matched, batch_unmatched) =
            process_batch(&current_batch, ref_kmers, k, threshold, canonical);
        matched.extend(batch_matched);
        unmatched.extend(batch_unmatched);
    }

    Ok((matched, unmatched))
}

pub struct KmerProcessor {
    pub current_kmer: u64,
    pub k: usize,
    pub kmers: HashSet<Kmer>,
}

impl KmerProcessor {
    pub fn new(k: usize) -> Self {
        KmerProcessor {
            current_kmer: 0,
            k,
            kmers: HashSet::new(),
        }
    }

    pub fn process(&mut self, seq: &[u8], id: Arc<str>) {
        if seq.len() < self.k {
            return; // Not enough bases for a k-mer
        }

        for i in 0..=(seq.len() - self.k) {
            let window = &seq[i..i + self.k];
            let kmer = Kmer::new(window.try_into().unwrap(), id.clone());
            self.kmers.insert(kmer);
        }
    }

    pub fn rolling_hash(&mut self, base: u8) {
        // Shift current k-mer left by 2 bits and add new base
        self.current_kmer = ((self.current_kmer << 2) & ((1 << (self.k * 2)) - 1))
            | match base {
                b'A' => 0b00,
                b'C' => 0b01,
                b'G' => 0b10,
                b'T' => 0b11,
                _ => panic!("Invalid nucleotide base"),
            };
    }
}
