use std::io::Read;

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

#[test]
fn test_kmer_fns() {
    let seq_vec = b"TGCTCAGATCATGTTTGTGTGG";
    let kmer = encode(seq_vec[0..21].try_into().unwrap());
    let kmer2 = encode(seq_vec[1..22].try_into().unwrap());

    println!("Encoded k-mer: {:042b}", kmer);
    println!("Encoded k-mer: {:042b}", kmer2);
    assert_ne!(kmer, kmer2);
    assert_eq!(kmer, 0b111001110100100011010011101111111011101110);

    let bit_cap = (1u64 << 21 * 2) - 1;
    let mut shifted_kmer = kmer;
    assert_eq!(shifted_kmer, kmer);

    shifted_kmer = ((kmer << 2) | encode(b"G")) & bit_cap;
    println!("Shifted k-mer: {:b}", shifted_kmer);

    assert_eq!(shifted_kmer, kmer2);
    assert_ne!(shifted_kmer, kmer);

    println!("Sequence vector: {:?}", seq_vec);

    let decoded_seq = [
        84, 71, 67, 84, 67, 65, 71, 65, 84, 67, 65, 84, 71, 84, 84, 84, 71, 84, 71, 84, 71, 71,
    ];

    assert_eq!(decode(kmer, 21), decoded_seq[0..21]);
    println!("Decoded k-mer 1: {:?}", decode(kmer, 21));

    assert_eq!(decode(kmer2, 21), decoded_seq[1..22]);
    println!("Decoded k-mer 2: {:?}", decode(kmer2, 21));

    assert_eq!(decode(shifted_kmer, 21), decoded_seq[1..22]);
    println!("Decoded shifter: {:?}", decode(shifted_kmer, 21));
}
