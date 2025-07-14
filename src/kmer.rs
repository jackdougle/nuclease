use std::collections::HashSet;

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
// pub fn extract_kmers(seq: &[u8], k: usize, canonical: bool) -> HashSet<Vec<u8>> {
//     let mut kmers = HashSet::new();

//     if seq.len() < k {
//         return kmers;
//     }

//     for i in 0..=(seq.len() - k) {
//         let window = &seq[i..i + k];

//         if window
//             .iter()
//             .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
//         {
//             continue; // skip kmers with ambiguous bases
//         }

//         let kmer = if canonical {
//             canonical_kmer(window)
//         } else {
//             window.to_ascii_uppercase()
//         };

//         kmers.insert(kmer);
//     }

//     kmers
// }
