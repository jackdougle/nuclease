use std::collections::HashSet;
use std::sync::Arc;
use ahash::AHashMap;
use std::hash::BuildHasherDefault;
use crate::kmer::*;

pub struct KmerProcessor {
    pub k: usize,
    pub current_kmers: HashSet<Kmer>,
    pub current_kmer: u64,
    pub threshold: usize,
    pub canonical: bool,
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: usize, canonical: bool) -> Self {
        KmerProcessor {
            k,
            current_kmers: HashSet::new(),
            current_kmer: 0,
            threshold,
            canonical,
        }
    }

    pub fn process_read(&mut self, read: &str) -> Result<(), String> {
        if read.len() < self.k {
            return Err("Read is shorter than k".to_string());
        }

        for i in 0..=read.len() - self.k {
            let kmer_seq = &read[i..i + self.k];
            let kmer = Kmer::new(kmer_seq.as_bytes().try_into().unwrap(), Arc::from(""));

            if self.canonical {
                let canonical_kmer = canonical_kmer(kmer_seq);
                self.current_kmers.insert(Kmer::new(canonical_kmer.as_bytes().try_into().unwrap(), Arc::from("")));
            } else {
                self.current_kmers.insert(kmer);
            }
        }

        Ok(())
    }
}