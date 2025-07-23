use std::process::exit;

use crate::kmer::*;
use rustc_hash::FxHashSet;

pub struct KmerProcessor {
    pub k: usize,
    pub threshold: u8,
    pub ref_kmers: FxHashSet<u64>,
    pub bit_cap: u64,
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: u8) -> Self {
        println!("KmerProcessor; k - {}, threshold - {}", k, threshold);
        KmerProcessor {
            k,
            threshold,
            ref_kmers: FxHashSet::default(),
            bit_cap: (1u64 << k * 2) - 1,
        }
    }

    pub fn process_ref(&mut self, ref_seq: Vec<u8>) {
        if ref_seq.len() < self.k {
            panic!("Read sequence is shorter than k");
        }

        let mut kmer: u64 = 0b00;

        for i in 0..=ref_seq.len() - self.k {
            if i == 0 {
                kmer = encode(&ref_seq[0..self.k]);
            } else {
                // Shift left by 2 bits and add the new base
                kmer = ((kmer << 2) | encode(&[ref_seq[i + self.k - 1]])) & self.bit_cap;
            }
            println!("Processed K-mer: {:?}", kmer);
            self.ref_kmers.insert(canonical_kmer(kmer, self.k));
        }
    }

    pub fn process_read(&self, read_seq: &str) -> bool {
        if read_seq.len() < self.k {
            println!(
                "Read sequence is shorter than k: {} < {}",
                read_seq.len(),
                self.k
            );
            return false;
        }

        let mut hits: u8 = 0;
        let mut kmer = 0b00;

        for i in 0..=read_seq.len() - self.k {
            if i == 0 {
                kmer = encode(&read_seq[0..self.k].as_bytes());
            } else {
                // Shift left by 2 bits and add the new base
                kmer = (kmer << 2) | encode(&[read_seq.as_bytes()[i + self.k - 1]]) & self.bit_cap;
            }
            if self.ref_kmers.contains(&canonical_kmer(kmer, self.k)) {
                hits += 1
            }
        }
        if hits >= self.threshold { true } else { false }
    }
}
