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
