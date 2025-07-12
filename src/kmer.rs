
pub fn reverse_complement(kmer: &[u8]) -> Vec<u8> {
    let mut rev_kmer = Vec::new();
    for &base in kmer.iter().rev() {
        rev_kmer.push(match base {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            _ => panic!("Invalid base: {}", base),
        });
    }
    rev_kmer
}