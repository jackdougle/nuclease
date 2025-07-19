use clap::Id;

use crate::Args;
use std::collections::HashSet;
use std::sync::Arc;

const K: usize = 21; // Default k-mer size, can be adjusted

pub struct IdKmer {
    pub id: Arc<str>,
    pub sequence: [u8; K], // Fixed size for k, input k-mer size, generally 21]
}

impl IdKmer {
    pub fn new(id: &str, sequence: [u8; K]) -> Self {
        IdKmer {
            id: Arc::from(id),
            sequence,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Kmer {
    pub encoded_seq: u64, // encoded using below definitions: 00011011 = ACGT
    pub id: Arc<str>,
}

impl Kmer {
    pub fn new(sequence: [u8; K], id: Arc<str>) -> Self {
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
            id: id,
        }
    }

    pub fn size_of(&self) -> usize {
        std::mem::size_of::<Self>() * 8 // size in bits
    }

    pub fn decoded(&self) -> [u8; K] {
        let mut seq = [0; K];
        let encoded = self.encoded_seq;
        for i in (0..K) {
            let base = match (encoded >> (i * 2)) & 0b11 {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _ => panic!("Non-DNA base!"),
            };
            seq[i] = base;
        }
        seq
    }
}

#[test]
fn test_kmer_struct() {
    let seq_vec = b"AGCTCAGATCATGTTTGTGTGG";
    let kmer = Kmer::new(seq_vec[0..K].try_into().unwrap(), Arc::from("SEQ_ID_1"));
    let kmer2 = Kmer::new(seq_vec[0..K].try_into().unwrap(), Arc::from("SEQ_ID_6"));

    println!("Encoded k-mer: {:042b}", kmer.encoded_seq);
    println!("Encoded k-mer: {:042b}", kmer2.encoded_seq);
    assert_eq!(
        kmer.encoded_seq,
        0b001001110100100011010011101111111011101110
    );

    assert_ne!(kmer, kmer2);
    assert_eq!(kmer.encoded_seq, kmer2.encoded_seq);
    println!(
        "k-mer 1: {}, k-mer 2: {}",
        kmer.encoded_seq, kmer2.encoded_seq
    );

    let decoded_seq: [u8; 21] = [
        71, 84, 71, 84, 71, 84, 84, 84, 71, 84, 65, 67, 84, 65, 71, 65, 67, 84, 67, 71, 65,
    ];

    println!("Decoded k-mer 1: {:?}", kmer.decoded());
    assert_eq!(kmer.decoded(), decoded_seq);

    let ex_str = b"AGCTCAGATCATGTTTGTGTG";

    // Convert [u8] to String (assuming ASCII ACGT bases)
    let decoded_str = String::from_utf8(decoded_seq.to_vec()).unwrap();
    println!("{}", decoded_str); // prints e.g. "AGTACGTGAC"
    assert_ne!(&kmer.decoded(), ex_str);
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
