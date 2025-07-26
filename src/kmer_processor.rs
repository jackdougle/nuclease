use needletail::sequence::canonical;
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

    pub fn process_ref(&mut self, ref_seq: &[u8]) {
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
            self.ref_kmers.insert(canonical_kmer(kmer, self.k));
        }
    }

    pub fn process_read(&self, read_seq: &[u8]) -> bool {
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
                kmer = encode(&read_seq[0..self.k]);
            } else {
                // Shift left by 2 bits and add the new base
                kmer = ((kmer << 2) | encode(&[read_seq[i + self.k - 1]])) & self.bit_cap;
            }
            //if self.ref_kmers.contains(&canonical_kmer(kmer, self.k)) {
            if self.ref_kmers.contains(&canonical_kmer(kmer, self.k)) {
                hits += 1;
                if hits >= self.threshold {
                    return true;
                }
            }
        }
        false
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

#[inline(always)]
pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut shift = (k - 1) * 2;
    for i in 0..k {
        let base = (kmer >> (i * 2)) & 0b11;
        let comp = base ^ 0b11;
        rc |= comp << shift;
        shift -= 2;
    }
    rc
}

#[inline(always)]
pub fn canonical_kmer(kmer: u64, k: usize) -> u64 {
    let rc = reverse_complement(kmer, k);
    if kmer < rc { kmer } else { rc }
}
