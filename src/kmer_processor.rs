use needletail::bitkmer::canonical;
use rustc_hash::FxHashSet;

pub struct KmerProcessor {
    pub k: usize,
    pub threshold: u8,
    pub ref_kmers: FxHashSet<u64>,
    pub bit_cap: u64,
}

impl KmerProcessor {
    pub fn new(k: usize, threshold: u8) -> Self {
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
                kmer = ((kmer << 2) | encode(&[ref_seq[i + self.k - 1]])) & self.bit_cap;
            }
            self.ref_kmers.insert(canonical((kmer, self.k as u8)).0.0);
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
                kmer = ((kmer << 2) | encode(&[read_seq[i + self.k - 1]])) & self.bit_cap;
            }

            if self
                .ref_kmers
                .contains(&canonical((kmer, self.k as u8)).0.0)
            {
                hits += 1;
                if hits >= self.threshold {
                    return true;
                }
            }
        }

        false
    }
}

#[inline(always)]
pub fn encode(seq: &[u8]) -> u64 {
    static BASE_TABLE: [u8; 85] = {
        let mut bases = [0u8; 85];
        bases[b'A' as usize] = 0b00;
        bases[b'C' as usize] = 0b01;
        bases[b'G' as usize] = 0b10;
        bases[b'T' as usize] = 0b11;
        bases
    };

    seq.iter().fold(0u64, |encoded, &base| {
        let val = BASE_TABLE[base as usize];
        debug_assert!(val <= 0b11, "Invalid base: {}", base as char);
        (encoded << 2) | val as u64
    })
}
