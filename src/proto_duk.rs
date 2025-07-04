use std::collections::HashMap;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{self, Read};
// use std::sync::{Arc, Mutex};
// use std::thread;
// use fastq::{Record, Parser};
// use std::io::{BufRead, BufReader, stdin};

fn read_file_to_string(filename: &str) -> Result<String, io::Error> {
    let mut file: File = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

fn load_seqs(path: &str) -> HashMap<String, String> {
    let mut seqs: HashMap<String, String> = HashMap::new();

    let seq_set = read_file_to_string(path)
        .expect(&format!("Error: Invalid reference file - {}", path));
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

fn rev_comp(seq: &str) -> String {
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

fn get_ref_kmers(ref_seqs: &HashMap<String, String>, k: usize, canonical_bool: bool) -> HashSet<String> {
    let mut ref_kmers: HashSet<String> = HashSet::new();

    for (_name, seq) in ref_seqs {
        for x in 0..=seq.len() - k {
            let temp = seq[x..x + k].to_string();

            if canonical_bool { 
                let rev = rev_comp(&temp);
                ref_kmers.insert(std::cmp::min(temp, rev));
            } else {
                ref_kmers.insert(temp);
            }

        }
    }

    ref_kmers
}

fn get_read_kmers(read_seqs: &HashMap<String, String>, k: usize) -> HashMap<String, Vec<String>> {
    let mut read_kmers: HashMap<String, Vec<String>> = HashMap::new();

    for (name, seq) in read_seqs {
        let expected_kmers = seq.len().saturating_sub(k) + 1;
        let kmers_vec = read_kmers.entry(name.clone()).or_insert_with(|| Vec::with_capacity(expected_kmers));
        
        for i in 0..=seq.len() - k {
            let kmer = &seq[i..i + k];
            
            kmers_vec.push(kmer.to_string());
        }
    }

    read_kmers
}

fn check_kmers(ref_kmers: &HashSet<String>, read_kmers: &HashMap<String, Vec<String>>, threshold: usize, canonical: bool) -> (HashMap<String, Vec<String>>, HashMap<String, Vec<String>>) {
    let mut matched: HashMap<String, Vec<String>> = HashMap::new();
    let mut unmatched: HashMap<String, Vec<String>> = HashMap::new();

    for (name, read_set) in read_kmers {
        let mut count = 0;
        let mut matched_kmers: Vec<String> = Vec::new();

        for read in read_set {
            let rev = rev_comp(read);
            let canonical = if canonical { std::cmp::min(read, &rev) } else { read };

            if ref_kmers.contains(canonical) { 
                count += 1;
                matched_kmers.push(canonical.clone());
            }
        }
        if count >= threshold {
            matched.entry(name.clone()).or_insert(Vec::new()).push(matched_kmers.join(", "));
        } else {
            unmatched.entry(name.clone()).or_insert(Vec::new()).push(read_set.join(", "));
        }
    }

    (matched, unmatched)
}

// TODO: Add threading infrastructure for parallel processing
pub fn main() {
    let args: Vec<String> = env::args().collect();
    let k = args[1].parse::<usize>().unwrap();
    let threshold = args[2].parse::<usize>().unwrap();
    let canonical = args[3].parse::<bool>().unwrap();

    let ref_filename = args[4].as_str();
    let ref_seqs = load_seqs(ref_filename);
    println!("ref_seqs: {:?}", ref_seqs);

    let read_filename = args[5].as_str();
    let read_seqs = load_seqs(read_filename);
    println!("read_seqs: {:?}", read_seqs);

    let ref_kmers = get_ref_kmers(&ref_seqs, k, canonical);
    let read_kmers = get_read_kmers(&read_seqs, k);

    let (matched, unmatched) = check_kmers(&ref_kmers, &read_kmers, threshold, canonical);

    println!("{}-mers found in reference database: {:?}", k, matched.len());
    println!("{}-mers not found in reference database: {:?}", k, unmatched.len());
}

