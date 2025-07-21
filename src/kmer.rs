use clap::Id;
use std::collections::HashSet;
use std::sync::Arc;

const K: usize = 21; // Default k-mer size, can be adjusted

#[test]
fn test_kmer_struct() {
    let seq_vec = b"AGCTCAGATCATGTTTGTGTGG";
    let kmer = encode(seq_vec[0..K].try_into().unwrap());
    let kmer2 = encode(seq_vec[0..K].try_into().unwrap());

    println!("Encoded k-mer: {:042b}", kmer);
    println!("Encoded k-mer: {:042b}", kmer2);
    assert_eq!(kmer, 0b001001110100100011010011101111111011101110);

    assert_ne!(kmer, kmer2);
    assert_eq!(kmer, kmer2);
    println!("k-mer 1: {}, k-mer 2: {}", kmer, kmer2);

    let decoded_seq: [u8; 21] = [
        71, 84, 71, 84, 71, 84, 84, 84, 71, 84, 65, 67, 84, 65, 71, 65, 67, 84, 67, 71, 65,
    ];

    println!("Decoded k-mer 1: {:?}", decode(kmer));
    assert_eq!(decode(kmer), decoded_seq);

    let ex_str = b"AGCTCAGATCATGTTTGTGTG";

    // Convert [u8] to String (assuming ASCII ACGT bases)
    let decoded_str = String::from_utf8(decoded_seq.to_vec()).unwrap();
    println!("{}", decoded_str); // prints e.g. "AGTACGTGAC"
    assert_ne!(&decode(kmer), ex_str);
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

pub fn decode(encoded: u64) -> [u8; K] {
    let mut seq = [0; K];
    for i in 0..K {
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

/// Return the reverse complement of a nucleotide sequence (A, C, G, T only).
pub fn reverse_complement_str(seq: &str) -> String {
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

pub fn reverse_complement(kmer: u64) -> u64 {
    let mut rc = 0u64;
    for i in 0..K {
        let base = (kmer >> (i * 2)) & 0b11;
        let comp = base ^ 0b11;
        rc |= comp << ((K - 1 - i) * 2);
    }
    rc
}

// Return the canonical form of a k-mer (lexicographically smallest of fwd and revcomp).
pub fn canonical_kmer(kmer: u64) -> u64 {
    let rc = reverse_complement(kmer);
    if kmer < rc { kmer } else { rc }
}
