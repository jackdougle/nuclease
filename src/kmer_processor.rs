use crate::kmer::*;
use ahash::AHashMap;
use needletail::bitkmer::canonical;
use std::collections::HashSet;
use std::hash::BuildHasherDefault;
use std::sync::Arc;

pub struct KmerProcessor {
    pub k: usize,
    pub threshold: u8,
    pub ref_kmers: HashSet<u64>,
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: u8) -> Self {
        KmerProcessor {
            k,
            threshold,
            ref_kmers: HashSet::new(),
        }
    }

    pub fn process_ref(&mut self, ref_seq: &str) {
        if ref_seq.len() < self.k {
            panic!("Read sequence is shorter than k");
        }

        for i in 0..=ref_seq.len() - self.k {}
    }

    pub fn process_read(&mut self, read_seq: &str) -> bool {
        if read_seq.len() < self.k {
            panic!("Reference sequence is shorter than k")
        }

        let mut count: u8 = 0;
        for i in 0..=read_seq.len() - self.k {
            let kmer = encode(&read_seq[i..i + self.k].as_bytes());
            let canonical = canonical_kmer(kmer);
            if self.ref_kmers.contains(&canonical) {
                count += 1
            }
        }
        if count >= self.threshold { true } else { false }
    }
}
