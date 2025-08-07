fn process_reads(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let processor = Arc::new(processor);

    let (tx, rx) = channel(); // for processed outputs
    let chunk_size = 1_000;
    let mut chunk = Vec::with_capacity(chunk_size);

    while let Some(record) = reader.next() {
        let record = record?;
        let id: Arc<str> = Arc::from(str::from_utf8(record.id())?);
        let seq = record.seq().into_owned();
        chunk.push((id, seq));

        if chunk.len() == chunk_size {
            let tx = tx.clone();
            let processor = processor.clone();
            let local_chunk = std::mem::take(&mut chunk);

            rayon::spawn(move || {
                let (matched, unmatched) = process_chunk(local_chunk, &processor);
                tx.send((matched, unmatched)).unwrap();
            });
        }
    }

    if !chunk.is_empty() {
        let tx = tx.clone();
        let processor = processor.clone();
        rayon::spawn(move || {
            let (matched, unmatched) = process_chunk(chunk, &processor);
            tx.send((matched, unmatched)).unwrap();
        });
    }

    drop(tx); // Close channel

    let mut matched_writer = BufWriter::new(File::create(matched_output)?);
    let mut unmatched_writer = BufWriter::new(File::create(unmatched_output)?);
    let mut mseq_count = 0;
    let mut mbase_count = 0;
    let mut useq_count = 0;
    let mut ubase_count = 0;

    for (matched, unmatched) in rx {
        mseq_count += matched.len() as u32;
        for (id, seq) in matched {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            mbase_count += seq.len() as u32;
        }

        useq_count += unmatched.len() as u32;
        for (id, seq) in unmatched {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            ubase_count += seq.len() as u32;
        }
    }

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_chunk(
    chunk: Vec<(Arc<str>, Vec<u8>)>,
    processor: &KmerProcessor,
) -> (Vec<(Arc<str>, Vec<u8>)>, Vec<(Arc<str>, Vec<u8>)>) {
    let mut matched = Vec::new();
    let mut unmatched = Vec::new();

    for (id, seq) in chunk {
        if processor.process_read(&seq) {
            matched.push((id, seq));
        } else {
            unmatched.push((id, seq));
        }
    }

    (matched, unmatched)
}
