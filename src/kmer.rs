use std::collections::HashSet;
use std::sync::Arc;

pub struct IdKmer {
    pub id: Arc<str>,
    pub sequence: [u8; 21], // Fixed size for k=21],
}

impl IdKmer {
    pub fn new(id: &str, sequence: [u8; 21]) -> Self {
        IdKmer {
            id: Arc::from(id),
            sequence,
        }
    }
}

pub struct Kmer {
    pub encoded_seq: u64, // encoded using below definitions: 00011011 = ACGT
}

impl Kmer {
    pub fn new(sequence: [u8; 21]) -> Self {
        Kmer {
            encoded_seq: sequence.iter().fold(0, |acc, &base| {
                (acc << 2)
                    | match base {
                        b'A' => 0b00,
                        b'C' => 0b01,
                        b'G' => 0b10,
                        b'T' => 0b11,
                        _ => panic!("Invalid nucleotide base"),
                    }
            }),
        }
    }

    pub fn size_of(&self) -> usize {
        std::mem::size_of::<Self>() * 8 // size in bits
    }
}

#[test]
fn test_kmer_creation() {
    let test_vec = b"ATGCTGTACGTAGCTCGATCA"; // 21 A's
    let kmer = Kmer::new(test_vec[..21].try_into().unwrap());
    println!("Encoded k-mer: {:042b}", kmer.encoded_seq);
    assert_eq!(kmer.size_of(), 64); // 64 bits for u64
}

/// Return the reverse complement of a nucleotide sequence (A, C, G, T only).
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|b| match b {
            'A' | 'a' => 'T',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'T' | 't' => 'A',
            _ => 'N', // default for non-ACGT
        })
        .collect()
}

/// Return the canonical form of a k-mer (lexicographically smallest of fwd and revcomp).
pub fn canonical_kmer(kmer: &str) -> String {
    let rc = reverse_complement(kmer);
    if kmer <= &rc {
        kmer.to_string()
    } else {
        rc.to_string()
    }
}

// Extract all k-mers from a sequence. Canonicalizes if requested.
// Skips over non-ACGT kmers (if they contain N or other symbols).
pub fn extract_kmers(seq: &[u8], k: usize, canonical: bool) -> HashSet<Vec<u8>> {
    let mut kmers = HashSet::new();

    if seq.len() < k {
        return kmers;
    }

    for i in 0..=(seq.len() - k) {
        let window = &seq[i..i + k];

        if window
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        {
            continue; // skip kmers with ambiguous bases
        }

        let kmer = if canonical {
            let temp = String::from_utf8(window.to_ascii_uppercase()).unwrap();
            canonical_kmer(&temp).as_bytes().to_vec()
        } else {
            window.to_ascii_uppercase()
        };

        kmers.insert(kmer);
    }

    kmers
}
