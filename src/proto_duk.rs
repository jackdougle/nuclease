use std::collections::HashMap;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{self, Read};
use std::process::exit;

fn read_file_to_string(filename: &str) -> Result<String, io::Error> {
    let mut file: File = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

fn load_seqs(path: &str) -> HashMap<String, String> {
    let mut seqs: HashMap<String, String> = HashMap::new();

    let seq_set = match read_file_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            println!("Error: Invalid reference file - {}", e);
            exit(1);
        }
    };
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

fn get_rev_comp(seq: &str) -> String {
    let mut rev: String = String::new();

    for c in seq.chars() {
        rev.push(match c {  
            'C' => 'G',
            'G' => 'C',
            'A' => 'T',
            'T' => 'A',
            _ => { println!("Invalid base detected: {}", c); c }
        });
    }

    rev.chars().rev().collect()
}

fn get_ref_kmers(ref_seqs: &HashMap<String, String>, k: usize) -> HashSet<String> {
    let mut ref_kmers: HashSet<String> = HashSet::new();

    for (_name, seq) in ref_seqs {
        let max = seq.len().saturating_sub(k.try_into().unwrap()) + 1;
        
        for x in 0..max {
            let end = x + k;

            let temp = seq[x..end].to_string();
            let rev = get_rev_comp(&temp);

            let canonical = std::cmp::min(temp, rev);

            ref_kmers.insert(canonical);
        }
    }

    ref_kmers
}

fn get_read_kmers(read_seqs: &HashMap<String, String>, k: usize) -> HashMap<String, Vec<String>> {
    let mut read_kmers: HashMap<String, Vec<String>> = HashMap::new();

    for (name, seq) in read_seqs {
        let max = seq.len().saturating_sub(k.try_into().unwrap()) + 1;

        for x in 0..max {
            let end = x + k;

            let temp = seq[x..end].to_string();

            read_kmers.entry(name.clone()).or_insert(Vec::new()).push(temp);
        }
    }

    read_kmers
}

fn check_kmers(ref_kmers: &HashSet<String>, read_kmers: &HashMap<String, Vec<String>>) -> HashMap<String, Vec<String>> {
    let mut sus_kmers: HashMap<String, Vec<String>> = HashMap::new();

    for (name, read_set) in read_kmers {
        for read in read_set {
            let rev = get_rev_comp(read);
            let canonical = std::cmp::min(read, &rev);

            if ref_kmers.contains(canonical) {
                sus_kmers.entry(name.clone()).or_insert(Vec::new()).push(read.to_string());
            }
        }
    }

    sus_kmers
}

pub fn main() {
    let args: Vec<String> = env::args().collect();
    let k = args[1].parse::<usize>().unwrap();

    let ref_filename = args[2].as_str();
    let ref_seqs = load_seqs(ref_filename);
    println!("ref_seqs: {:?}", ref_seqs);

    let read_filename = args[3].as_str();
    let read_seqs = load_seqs(read_filename);
    println!("read_seqs: {:?}", read_seqs);

    let ref_kmers = get_ref_kmers(&ref_seqs, k);
    let read_kmers = get_read_kmers(&read_seqs, k);

    let sus_kmers = check_kmers(&ref_kmers, &read_kmers);

    println!("{}-mers found in reference database: {:?}", k, sus_kmers);
}