use crate::kmer_processor::encode;
use std::time::Instant;
use wide::u8x16;

pub fn encode_simd(seq: &[u8]) -> u64 {
    assert!(seq.len() <= 32);

    if seq.len() <= 16 {
        // Pad with 'A' (0) to 16 if needed
        let mut padded = [b'A'; 16];
        padded[..seq.len()].copy_from_slice(seq);
        let packed = encode_16(padded);
        // Shift right to discard padded bases at the right
        packed >> (2 * (16 - seq.len()))
    } else {
        // length between 17 and 32
        let mut arr: [u8; 16] = [b'A'; 16];
        arr[..16].copy_from_slice(&seq[..16]);
        let high = encode_16(arr);
        let low = encode_simd(&seq[16..]); // recursive call for rest (<=16)
        (high << (2 * (seq.len() - 16))) | low
    }
}

pub fn encode_16(chunk: [u8; 16]) -> u64 {
    // Load 16 bytes into u8x16
    let v = u8x16::new(chunk);

    // Create masks for bases
    let is_c = v.cmp_eq(u8x16::splat(b'C'));
    let is_g = v.cmp_eq(u8x16::splat(b'G'));
    let is_t = v.cmp_eq(u8x16::splat(b'T'));

    // Map bases to 2-bit codes
    // code = (is_c * 1) | (is_g * 2) | (is_t * 3), 'A' is 0
    // wide::Mask<T> can be cast to u8x16 with .to_int()
    let codes = (is_c & u8x16::splat(1)) | (is_g & u8x16::splat(2)) | (is_t & u8x16::splat(3));

    // Extract lanes as array for bit packing
    let lanes = codes.to_array();

    // Pack 16 2-bit codes into u64 (MSB = first base)
    let mut packed: u64 = 0;
    for &code in &lanes {
        packed = (packed << 2) | (code as u64);
    }

    packed
}

#[test]
fn simd_test() {
    let timer = Instant::now();
    let seq: [u8; 16] = *b"ACGTACGTACGTACGT";
    let packed = encode_16(seq);
    println!("{:042b}", packed);
    println!("Time:\t{}", timer.elapsed().as_nanos());
}

#[test]
fn sisd_test() {
    let timer = Instant::now();
    let seq = b"ACGTACGTACGTACGT";
    let packed = encode(seq);
    println!("{:042b}", packed);
    println!("Time:\t{}", timer.elapsed().as_nanos());
}

#[test]
fn compare() {
    let sisd_out = encode(b"ACGTACGTACGTACGT");
    let simd_out = encode_simd(b"ACGTACGTACGTACGT");
    assert_eq!(sisd_out, simd_out);
}
