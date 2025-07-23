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

pub fn encode(sequence: &[u8]) -> u64 {
    sequence.iter().fold(0, |acc, &base| {
        (acc << 2)
            | match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => panic!("Invalid nucleotide base"),
            }
    })
}

pub fn decode(encoded: u64, k: usize) -> Vec<u8> {
    let mut seq = Vec::new();
    for i in 0..k {
        let base = match (encoded >> (i * 2)) & 0b11 {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => panic!("Non-DNA base!"),
        };
        seq.push(base);
    }
    seq.into_iter().rev().collect()
}

pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for i in 0..k {
        let base = (kmer >> (i * 2)) & 0b11;
        let comp = base ^ 0b11;
        rc |= comp << ((k - 1 - i) * 2);
    }
    rc
}

// Return the canonical form of a k-mer (lexicographically smallest of fwd and revcomp).
pub fn canonical_kmer(kmer: u64, k: usize) -> u64 {
    let rc = reverse_complement(kmer, k);
    if kmer < rc { kmer } else { rc }
}
