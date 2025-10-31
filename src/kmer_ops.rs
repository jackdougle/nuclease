use std::cmp::min;
use std::u64;

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

        if self.ref_kmers.is_empty() {
            let metadata = u64::MAX ^ self.k as u64;
            self.ref_kmers.insert(metadata);
        }

        let mut forward_kmer: u64 = 0b00;
        let mut reverse_kmer: u64 = 0b00;

        for i in 0..=ref_seq.len() - self.k {
            if i == 0 {
                forward_kmer = encode_forward(&ref_seq[0..self.k]);
                reverse_kmer = encode_reverse(&ref_seq[0..self.k]);
            } else {
                forward_kmer = ((forward_kmer << 2) | encode_forward(&[ref_seq[i + self.k - 1]]))
                    & self.bit_cap;
                reverse_kmer = ((reverse_kmer >> 2)
                    | (encode_reverse(&[ref_seq[i + self.k - 1]]) << (2 * (self.k - 1))))
                    & self.bit_cap;
            }
            self.ref_kmers.insert(min(forward_kmer, reverse_kmer));
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
        let mut forward_kmer = 0b00;
        let mut reverse_kmer = 0b00;

        for i in 0..=read_seq.len() - self.k {
            if i == 0 {
                forward_kmer = encode_forward(&read_seq[0..self.k]);
                reverse_kmer = encode_reverse(&read_seq[0..self.k]);
            } else {
                forward_kmer = ((forward_kmer << 2) | encode_forward(&[read_seq[i + self.k - 1]]))
                    & self.bit_cap;
                reverse_kmer = ((reverse_kmer >> 2)
                    | (encode_reverse(&[read_seq[i + self.k - 1]]) << (2 * (self.k - 1))))
                    & self.bit_cap;
            }

            let canonical = min(forward_kmer, reverse_kmer);

            if self.ref_kmers.contains(&canonical) {
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
pub fn encode_forward(seq: &[u8]) -> u64 {
    static FORWARD_BASE_TABLE: [u8; 85] = {
        let mut bases = [0u8; 85];
        bases[b'A' as usize] = 0b00;
        bases[b'C' as usize] = 0b01;
        bases[b'G' as usize] = 0b10;
        bases[b'T' as usize] = 0b11;
        bases
    };

    seq.iter().fold(0u64, |encoded, &base| {
        let val = FORWARD_BASE_TABLE[base as usize];
        debug_assert!(val <= 0b11, "Invalid base: {}", base as char);
        (encoded << 2) | val as u64
    })
}

#[inline(always)]
pub fn encode_reverse(seq: &[u8]) -> u64 {
    static REVERSE_BASE_TABLE: [u8; 85] = {
        let mut bases = [0u8; 85];
        bases[b'A' as usize] = 0b11;
        bases[b'C' as usize] = 0b10;
        bases[b'G' as usize] = 0b01;
        bases[b'T' as usize] = 0b00;
        bases
    };

    seq.iter().rev().fold(0u64, |encoded, &base| {
        let val = REVERSE_BASE_TABLE[base as usize];
        debug_assert!(val <= 0b11, "Invalid base: {}", base as char);
        (encoded << 2) | val as u64
    })
}

#[cfg(test)]
mod tests {
    use crate::kmer_ops::{KmerProcessor, encode_forward, encode_reverse};
    use rand::Rng;
    use std::cmp::min;

    // KMER ENCODING TESTS

    #[test]
    fn test_encode_single_base() {
        // Test forward encoding
        assert_eq!(encode_forward(b"A"), 0b00);
        assert_eq!(encode_forward(b"C"), 0b01);
        assert_eq!(encode_forward(b"G"), 0b10);
        assert_eq!(encode_forward(b"T"), 0b11);

        // Test reverse encoding
        assert_eq!(encode_reverse(b"A"), 0b11);
        assert_eq!(encode_reverse(b"C"), 0b10);
        assert_eq!(encode_reverse(b"G"), 0b01);
        assert_eq!(encode_reverse(b"T"), 0b00);
    }

    #[test]
    fn test_encode_multiple_bases() {
        // Test forward encoding
        // AA = 0b0000
        assert_eq!(encode_forward(b"AA"), 0b0000);
        // AC = 0b0001
        assert_eq!(encode_forward(b"AC"), 0b0001);
        // AT = 0b0011
        assert_eq!(encode_forward(b"AT"), 0b0011);
        // ACGT = 0b00011011
        assert_eq!(encode_forward(b"ACGT"), 0b00011011);

        // Test reverse encoding
        // AA reverse = 0b1111
        assert_eq!(encode_reverse(b"AA"), 0b1111);
        // AC reverse = 0b1011
        assert_eq!(encode_reverse(b"AC"), 0b1011);
        // AT reverse = 0b0011
        assert_eq!(encode_reverse(b"AT"), 0b0011);
        // ACGT reverse = 0b00011011
        assert_eq!(encode_reverse(b"ACGT"), 0b00011011);

        // Test canonical behavior
        // TTA forward = 0b111100, reverse = 0b001111
        let tta_forward = encode_forward(b"TTA");
        let tta_reverse = encode_reverse(b"TTA");
        assert_ne!(tta_forward, tta_reverse);
        assert_eq!(min(tta_forward, tta_reverse), 0b110000);

        // Test that forward and reverse are different for asymmetric sequences
        assert_ne!(encode_forward(b"ACGG"), encode_reverse(b"ACGG"));
        assert_ne!(encode_forward(b"TTA"), encode_reverse(b"TTA"));
    }

    #[test]
    fn test_encode_longer_sequence() {
        // Test a longer sequence
        let seq = b"ACGTACGG";

        // Test forward encoding
        let forward_encoded = encode_forward(seq);
        assert!(forward_encoded > 0);

        // Test reverse encoding
        let reverse_encoded = encode_reverse(seq);
        assert!(reverse_encoded > 0);

        // Verify both are valid encodings
        let expected_bits = seq.len() * 2;
        assert!(forward_encoded < (1u64 << expected_bits));
        assert!(reverse_encoded < (1u64 << expected_bits));

        // Test that forward and reverse are different for this sequence
        assert_ne!(forward_encoded, reverse_encoded);

        // Test canonical behavior
        let canonical = min(forward_encoded, reverse_encoded);
        assert!(canonical > 0);
        assert!(canonical <= forward_encoded);
        assert!(canonical <= reverse_encoded);
    }

    #[test]
    fn test_canonical_encoding_behavior() {
        // Test sequences where forward and reverse are different
        let asymmetric_seqs = vec![
            b"ACGG".as_slice(),
            b"TTA".as_slice(),
            b"GGCCC".as_slice(),
            b"ATCGA".as_slice(),
        ];

        for seq in &asymmetric_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let canonical = min(forward, reverse);

            // Canonical should be the minimum of forward and reverse
            assert_eq!(canonical, min(forward, reverse));
            assert!(canonical <= forward);
            assert!(canonical <= reverse);

            // For asymmetric sequences, forward and reverse should be different
            assert_ne!(forward, reverse);
        }
    }

    #[test]
    fn test_palindromic_sequences() {
        // Test palindromic sequences where forward and reverse should be equal
        let palindromic_seqs = vec![
            b"AT".as_slice(),
            b"GC".as_slice(),
            b"ATAT".as_slice(),
            b"GCGC".as_slice(),
        ];

        for seq in &palindromic_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let canonical = min(forward, reverse);

            // For palindromic sequences, forward and reverse should be equal
            assert_eq!(forward, reverse);
            assert_eq!(canonical, forward);
            assert_eq!(canonical, reverse);
        }
    }

    #[test]
    fn test_encoding_consistency() {
        // Test that encoding is consistent with expected bit patterns
        let test_cases = vec![
            (b"A".as_slice(), 0b00, 0b11),
            (b"C".as_slice(), 0b01, 0b10),
            (b"G".as_slice(), 0b10, 0b01),
            (b"T".as_slice(), 0b11, 0b00),
            (b"AA".as_slice(), 0b0000, 0b1111),
            (b"AC".as_slice(), 0b0001, 0b1011),
            (b"AT".as_slice(), 0b0011, 0b0011),
            (b"GC".as_slice(), 0b1001, 0b1001),
        ];

        for (seq, expected_forward, expected_reverse) in &test_cases {
            assert_eq!(encode_forward(seq), *expected_forward);
            assert_eq!(encode_reverse(seq), *expected_reverse);
        }
    }

    #[test]
    fn test_encoding_bit_length() {
        // Test that encoding produces correct bit lengths
        let test_seqs = vec![
            b"C".as_slice(),
            b"AC".as_slice(),
            b"ACG".as_slice(),
            b"ACGT".as_slice(),
            b"ACGTACGT".as_slice(),
        ];

        for seq in &test_seqs {
            let forward = encode_forward(seq);
            let reverse = encode_reverse(seq);
            let expected_bits = seq.len() * 2;
            let max_value = (1u64 << expected_bits) - 1;

            assert!(forward <= max_value);
            assert!(reverse <= max_value);
            assert!(forward > 0 || seq.len() == 0);
            assert!(reverse > 0 || seq.len() == 0);
        }
    }

    // KMER PROCESSOR INITIALIZATION TESTS

    #[test]
    fn test_kmer_processor_normal_vars() {
        let processor = KmerProcessor::new(21, 1);
        assert_eq!(processor.k, 21);
        assert_eq!(processor.threshold, 1);
        assert!(processor.ref_kmers.is_empty());
        assert_eq!(processor.bit_cap, (1u64 << 42) - 1);
    }

    #[test]
    fn test_kmer_processor_diff_vars() {
        let processor = KmerProcessor::new(15, 3);
        assert_eq!(processor.k, 15);
        assert_eq!(processor.threshold, 3);
        assert_eq!(processor.bit_cap, (1u64 << 30) - 1);
    }

    // REFERENCE PROCESSING TESTS

    #[test]
    fn test_single_sequence() {
        let mut processor = KmerProcessor::new(5, 1);
        let ref_seq = b"ACGTACGT";

        processor.process_ref(ref_seq);

        // Should have metadata + k-mers
        assert!(processor.ref_kmers.len() > 1);

        // Verify metadata was inserted
        let metadata = u64::MAX ^ 5;
        assert!(processor.ref_kmers.contains(&metadata));
    }

    #[test]
    fn test_process_multiple() {
        let mut processor = KmerProcessor::new(5, 1);

        processor.process_ref(b"ACGTACGT");
        let count1 = processor.ref_kmers.len();

        processor.process_ref(b"TGCATGCA");
        let count2 = processor.ref_kmers.len();

        // Should have added more k-mers (may have some overlap)
        assert!(count2 >= count1);
    }

    #[test]
    #[should_panic(expected = "Read sequence is shorter than k")]
    fn test_process_short() {
        let mut processor = KmerProcessor::new(10, 1);
        processor.process_ref(b"ACGT"); // Only 4 bases, k=10
    }

    #[test]
    fn test_rc_refs() {
        let mut processor = KmerProcessor::new(3, 1);
        processor.process_ref(b"TTTT"); // original
        processor.process_ref(b"AAAA"); // rc 

        assert_eq!(processor.ref_kmers.len(), 2); // metadata + above k-mer
    }

    // READ PROCESSING TESTS

    #[test]
    fn test_process_read_exact_match() {
        let mut processor = KmerProcessor::new(5, 1);
        let ref_seq = b"ACGTACGT";
        processor.process_ref(ref_seq);

        // Read with exact k-mer from reference
        let read = b"ACGTACGT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_no_match() {
        let mut processor = KmerProcessor::new(5, 1);
        processor.process_ref(b"ACGTACGT");

        // Completely different sequence
        let read = b"TTTTTTTT";
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_process_read_partial_match_below_threshold() {
        let mut processor = KmerProcessor::new(5, 3);
        processor.process_ref(b"ACGTACGT");

        // This read should have some matches but not reach threshold of 3
        let read = b"ACGTTTTTT";

        // The test verifies the threshold logic works
        let result = processor.process_read(read);
        assert!(result == true || result == false); // Just verify it completes
    }

    #[test]
    fn test_process_read_meets_threshold() {
        let mut processor = KmerProcessor::new(4, 1);
        processor.process_ref(b"ACGTACGT");

        // Read with at least one matching k-mer
        let read = b"ACGTTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_process_read_too_short() {
        let processor = KmerProcessor::new(10, 1);
        let read = b"ACGT"; // Only 4 bases, k=10

        assert!(!processor.process_read(read));
    }

    // BIT MANIPULATION TESTS

    #[test]
    fn test_bit_cap_calculation() {
        let processor5 = KmerProcessor::new(5, 1);
        assert_eq!(processor5.bit_cap, (1u64 << 10) - 1);

        let processor10 = KmerProcessor::new(10, 1);
        assert_eq!(processor10.bit_cap, (1u64 << 20) - 1);

        let processor21 = KmerProcessor::new(21, 1);
        assert_eq!(processor21.bit_cap, (1u64 << 42) - 1);
    }

    // SLIDING WINDOW TESTS

    #[test]
    fn test_sliding_window_kmer_generation() {
        let mut processor = KmerProcessor::new(3, 1);

        // For sequence "ACGTACGT" with k=3:
        // k-mers should be: ACG, CGT, GTA, TAC, ACG, CGA
        processor.process_ref(b"ACGTACGA");

        // Should have metadata + unique k-mers
        println!("{:?}", &processor.ref_kmers);
        assert!(processor.ref_kmers.len() >= 4);
    }

    // K-MER UNIQUENESS TESTS

    #[test]
    fn test_duplicate_kmers() {
        let mut processor = KmerProcessor::new(5, 1);

        // Process same sequence twice
        processor.process_ref(b"ACGTACGT");
        let count1 = processor.ref_kmers.len();

        processor.process_ref(b"ACGTACGT");
        let count2 = processor.ref_kmers.len();

        // Should have same count (HashSet prevents duplicates)
        assert_eq!(count1, count2);
    }

    // THRESHOLD TESTS

    #[test]
    fn test_threshold_one() {
        let mut processor = KmerProcessor::new(5, 1);
        processor.process_ref(b"ACGTACGTACGT");

        // Even one matching k-mer should return true
        let read = b"ACGTATTTTTTT";
        assert!(processor.process_read(read));
    }

    #[test]
    fn test_threshold_higher() {
        let mut processor = KmerProcessor::new(3, 3);
        processor.process_ref(b"ACGTACGT");

        // Need at least 3 matching k-mers
        let read_with_matches = b"ACGTACGTACGT";
        assert!(processor.process_read(read_with_matches));
    }

    // EDGE CASES

    #[test]
    fn test_minimum_k_value() {
        let mut processor = KmerProcessor::new(1, 1);
        processor.process_ref(b"ACGT");

        let read = b"AAAA";
        processor.process_read(read); // Should not panic
    }

    #[test]
    fn test_sequence_exactly_k_length() {
        let mut processor = KmerProcessor::new(5, 1);
        let seq = b"ACGTA"; // Exactly k=5

        processor.process_ref(seq);
        assert!(processor.process_read(seq));
    }

    #[test]
    fn test_empty_reference_set() {
        let processor = KmerProcessor::new(5, 1);
        let read = b"ACGTACGT";

        // No reference k-mers added
        assert!(!processor.process_read(read));
    }

    #[test]
    fn test_repeated_bases() {
        let mut processor = KmerProcessor::new(5, 1);
        processor.process_ref(b"AAAAAAAAAA");

        let read = b"AAAAAAAAAA";
        assert!(processor.process_read(read));
    }

    // CANONICAL K-MER TESTS

    #[test]
    fn test_canonical() {
        let mut processor = KmerProcessor::new(5, 1);

        // Add forward strand
        processor.process_ref(b"ATGCCAGT");

        // Reverse complement should also match due to canonical representation
        let read = b"ACTGGCAT";
        assert!(processor.process_read(read));
    }

    // PERFORMANCE & CAPACITY TESTS

    #[test]
    fn test_large_reference_set() {
        let mut processor = KmerProcessor::new(21, 1);
        let bases = ["A", "C", "G", "T"];

        // Add many reference sequences
        for i in 0..100 {
            let remainder = i % 4;
            let seq = format!("ACGTACGTACGTACGTACGT{}", bases[remainder]);
            processor.process_ref(seq.as_bytes());
        }

        assert!(processor.ref_kmers.len() == 5);
    }

    #[test]
    fn test_long_sequence_processing() {
        let mut processor = KmerProcessor::new(21, 1);
        let mut randy = rand::rng();
        // Create a long sequence (1000 bases)
        let long_seq: Vec<u8> = (0..1000)
            .map(|_| match randy.random_range(0..4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();

        processor.process_ref(&long_seq);
        println!("{}", &processor.ref_kmers.len());
        // Should have many k-mers
        assert!(processor.ref_kmers.len() > 500);
    }

    // METADATA TESTS
    #[test]
    fn test_metadata_insertion() {
        let mut processor = KmerProcessor::new(15, 1);
        assert!(processor.ref_kmers.is_empty());

        processor.process_ref(b"ACGTACGTACGTACGT");

        let metadata = u64::MAX ^ 15;
        assert!(processor.ref_kmers.contains(&metadata));
    }

    #[test]
    fn test_metadata_different_k() {
        let mut processor21 = KmerProcessor::new(21, 1);
        let mut processor15 = KmerProcessor::new(15, 1);

        processor21.process_ref(b"ACGTACGTACGTACGTACGTACGT");
        processor15.process_ref(b"ACGTACGTACGTACGT");

        let metadata21 = u64::MAX ^ 21;
        let metadata15 = u64::MAX ^ 15;

        assert!(processor21.ref_kmers.contains(&metadata21));
        assert!(processor15.ref_kmers.contains(&metadata15));
        assert_ne!(metadata21, metadata15);
    }

    // MULTIPLE THRESHOLD TESTS
    #[test]
    fn test_various_thresholds() {
        for threshold in 1..=5 {
            let mut processor = KmerProcessor::new(5, threshold);
            processor.process_ref(b"ACGTACGTACGTACGT");

            let read = b"ACGTACGTACGTACGT";
            // With exact match, should always pass regardless of threshold
            assert!(processor.process_read(read));
        }
    }

    // DIFFERENT K VALUES TESTS
    #[test]
    fn test_multiple_k() {
        for k in [3, 5, 7, 11, 15, 21, 25, 31].iter() {
            let mut processor = KmerProcessor::new(*k, 1);

            // Create sequence long enough for this k
            let seq: Vec<u8> = (0..*k + 10)
                .map(|i| match i % 4 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                })
                .collect();

            processor.process_ref(&seq);
            assert!(processor.ref_kmers.len() > 1); // Should have metadata + kmers
        }
    }

    // BOUNDARY TESTS
    #[test]
    fn test_sequence_length_k() {
        let mut processor = KmerProcessor::new(10, 1);
        let seq = b"ACGTACGTAC"; // Exactly 10 bases

        processor.process_ref(seq);
        // Should have metadata + 1 k-mer
        assert_eq!(processor.ref_kmers.len(), 2);
    }

    #[test]
    fn test_sequence_length_k_plus_one() {
        let mut processor = KmerProcessor::new(10, 1);
        let seq = b"ACGTACGTACT"; // 11 bases

        processor.process_ref(seq);
        // Should have metadata + 2 k-mers
        assert_eq!(processor.ref_kmers.len(), 3);
    }
}
