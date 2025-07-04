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
