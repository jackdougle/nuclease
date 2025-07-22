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
    seq
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
fn test_kmer_struct() {
    let seq_vec = b"AGCTCAGATCATGTTTGTGTGG";
    let kmer = encode(seq_vec[0..21].try_into().unwrap());
    let kmer2 = encode(seq_vec[0..21].try_into().unwrap());

    println!("Encoded k-mer: {:042b}", kmer);
    println!("Encoded k-mer: {:042b}", kmer2);
    assert_eq!(kmer, 0b001001110100100011010011101111111011101110);

    assert_ne!(kmer, kmer2);
    assert_eq!(kmer, kmer2);
    println!("k-mer 1: {}, k-mer 2: {}", kmer, kmer2);

    let decoded_seq: [u8; 21] = [
        71, 84, 71, 84, 71, 84, 84, 84, 71, 84, 65, 67, 84, 65, 71, 65, 67, 84, 67, 71, 65,
    ];

    println!("Decoded k-mer 1: {:?}", decode(kmer, 21));
    assert_eq!(decode(kmer, 21), decoded_seq);

    let ex_str = b"AGCTCAGATCATGTTTGTGTG";

    // Convert [u8] to String (assuming ASCII ACGT bases)
    let decoded_str = String::from_utf8(decoded_seq.to_vec()).unwrap();
    println!("{}", decoded_str); // prints e.g. "AGTACGTGAC"
    assert_ne!(&decode(kmer, 21), ex_str);
}
