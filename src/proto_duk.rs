use std::collections::HashMap;
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

fn load_seqs(path: &str, k: usize, set_type: bool) -> HashMap<String, String> {
    let mut seqs: HashMap<String, String> = HashMap::new();

    println!("Loading seqs from {}", set_type);

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

fn get_kmers(seqs: &HashMap<String, String>, k: usize, canonical: bool) -> HashMap<String, String> {
    let mut kmers: HashMap<String, String> = HashMap::new();

    for (name, seq) in seqs.iter() {
        let max: usize = seq.len().saturating_sub(k.try_into().unwrap()) + 1;

        for x in 0..max {
            println!("{}", x);
            let end: usize = x.checked_add(k.try_into().unwrap()).unwrap_or(x);
            let temp_seq = &seq[x..end];
            kmers.insert(name.to_string(), temp_seq.to_string());
            if canonical {
                let rev = get_complement(temp_seq);
                kmers.insert(name.to_string(), rev.clone());
                if x == 50 {
                    println!("Here!");
                    println!("{}", rev);
                    println!("{}", temp_seq);
                }
            }
        }
    }

    kmers
}

// fn check_kmers<'a>(k: i32, tests: HashMap<&'a str, &'a str>) -> HashMap<&'a str, &'a str> {
//     let mut sus: HashMap<&'a str, &'a str> = HashMap::new();

//     for (name, seq) in tests {
//         let max: usize = seq.len().saturating_sub(k.try_into().unwrap()) + 1;
//         println!("{}", seq);

//         for x in 0..max {
//             let end: usize = x.checked_add(k.try_into().unwrap()).unwrap_or(x);
//             let temp_seq = &seq[x..end];
//             println!("start: {}, end: {}, {}-mer: {:?}", x, end, k, temp_seq);
//         }
//     }

//     sus
// }

fn get_complement(seq: &str) -> String {
    let mut rev: String = String::new();

    for c in seq.chars() {
        if c == 'C' {
            rev.push('G');
        } else if c == 'G' {
            rev.push('C');
        } else if c == 'A' {
            rev.push('T');
        } else if c == 'T' {
            rev.push('A');
        } else {
            println!("Invalid base detected");
        }
    }

    rev
}

pub fn main() {
    let args: Vec<String> = env::args().collect();
    let k = args[1].parse::<usize>().unwrap();

    let ref_seqs = load_seqs("input/refs.fasta", k, true);
    println!("ref_seqs: {:?}", ref_seqs);
    let test_seqs = load_seqs("input/tests.fastq", k, true);
    println!("test_seqs: {:?}", test_seqs);

    let ref_kmers = get_kmers(&ref_seqs, k, true);
    let test_kmers = get_kmers(&test_seqs, k, false);

    // TODO: make function to get kmers from seqs
    // TODO: make function to get complements of kmers
    // TODO: make function to check sample kmers against refs

    //let sus: HashMap<&str, &str> = check_kmers(k, samples);

    //let full_refs: HashMap<String, String> = get_complements(&refs);
}
