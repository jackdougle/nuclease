pub fn canonical_kmer(seq: &str) -> String {
    let rc: String = seq
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect();

    if rc.as_str() < seq {
        rc
    } else {
        seq.to_string()
    }
}
