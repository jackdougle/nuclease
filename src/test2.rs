fn process_reads_old(
    reads_path: &str,
    processor: KmerProcessor,
    matched_path: &str,
    unmatched_path: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let processor = Arc::new(processor);

    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_path)?)));
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_path)?)));

    let mut chunks = Vec::new();
    let mut chunk = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;
        let id = Arc::from(str::from_utf8(record.id())?);
        let seq = record.seq().into_owned();
        chunk.push((id, seq));

        if chunk.len() == 10_000 {
            chunks.push(chunk);
            chunk = Vec::new();
        }
    }

    if !chunk.is_empty() {
        chunks.push(chunk);
    }

    // Process all chunks in parallel and collect their results
    let results: Vec<(u32, u32, u32, u32)> = chunks
        .into_par_iter()
        .map(|chunk| {
            process_chunk(
                chunk,
                processor.clone(),
                matched_writer.clone(),
                unmatched_writer.clone(),
            )
        })
        .collect::<Result<_, _>>()?;

    // Reduce the results (sum them)
    let (mseq_count, mbase_count, useq_count, ubase_count) =
        results
            .into_iter()
            .fold((0u32, 0u32, 0u32, 0u32), |total, local| {
                (
                    total.0 + local.0,
                    total.1 + local.1,
                    total.2 + local.2,
                    total.3 + local.3,
                )
            });

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_chunk(
    reads_chunk: Vec<(Arc<str>, Vec<u8>)>,
    processor: Arc<KmerProcessor>,
    matched_writer: Arc<Mutex<BufWriter<File>>>,
    unmatched_writer: Arc<Mutex<BufWriter<File>>>,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let results: Vec<(bool, Arc<str>, Vec<u8>)> = reads_chunk
        .into_par_iter()
        .map(|(id, seq)| {
            let is_match = processor.process_read(&seq);
            (is_match, id, seq)
        })
        .collect();

    let mut matched_writer = matched_writer.lock().unwrap();
    let mut unmatched_writer = unmatched_writer.lock().unwrap();

    let (mut mseq_count, mut mbase_count) = (0u32, 0u32);
    let (mut useq_count, mut ubase_count) = (0u32, 0u32);

    for (is_match, id, seq) in results {
        if is_match {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            mseq_count += 1;
            mbase_count += seq.len() as u32;
        } else {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            useq_count += 1;
            ubase_count += seq.len() as u32;
        }
    }

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_reads_crossbeam(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let processor = Arc::new(processor);

    let (tx, rx) = channel(); // for processed outputs
    let chunk_size = 10_000;
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
                let (matched, unmatched) = process_chunk_crossbeam(local_chunk, &processor);
                tx.send((matched, unmatched)).unwrap();
            });
        }
    }

    if !chunk.is_empty() {
        let tx = tx.clone();
        let processor = processor.clone();
        rayon::spawn(move || {
            let (matched, unmatched) = process_chunk_crossbeam(chunk, &processor);
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
        for (id, seq) in matched {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            mseq_count += 1;
            mbase_count += seq.len() as u32;
        }
        for (id, seq) in unmatched {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            useq_count += 1;
            ubase_count += seq.len() as u32;
        }
    }

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_chunk_crossbeam(
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

// Bad print statements, first Rayon architecture
fn get_read_chunks(
    reads_path: &str,
    processor: KmerProcessor,
    matched_output: &str,
    unmatched_output: &str,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let mut chunk = Vec::new();

    let processor = Arc::new(processor);

    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_output)?)));
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_output)?)));

    let matched_count = Arc::new(Mutex::new(0u32));
    let unmatched_count = Arc::new(Mutex::new(0u32));

    let mut chunks: Vec<Vec<(String, Vec<u8>)>> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;
        let id = String::from_utf8(Vec::from(record.id()))?;
        let seq = record.seq().into_owned();

        chunk.push((id, seq));

        if chunk.len() == 10_000 {
            chunks.push(chunk);
            chunk = Vec::new();
        }
    }

    if !chunk.is_empty() {
        chunks.push(chunk);
    }

    chunks.into_par_iter().try_for_each(|chunk| {
        process_reads_parallel(
            chunk,
            processor.clone(),
            matched_writer.clone(),
            unmatched_writer.clone(),
            matched_count.clone(),
            unmatched_count.clone(),
        )
    })?;

    println!(
        "Wrote {} matched reads to {}!",
        matched_count.lock().unwrap(),
        matched_output
    );
    println!(
        "Wrote {} unmatched reads to {}!",
        unmatched_count.lock().unwrap(),
        unmatched_output
    );

    Ok(())
}

fn process_reads_parallel(
    reads_chunk: Vec<(String, Vec<u8>)>,
    processor: Arc<KmerProcessor>,
    matched_writer: Arc<Mutex<BufWriter<File>>>,
    unmatched_writer: Arc<Mutex<BufWriter<File>>>,
    matched_count: Arc<Mutex<u32>>,
    unmatched_count: Arc<Mutex<u32>>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let results: Vec<(bool, String, Vec<u8>)> = reads_chunk
        .into_par_iter()
        .map(|(id, seq)| {
            let is_match = processor.process_read(&seq);
            (is_match, id, seq)
        })
        .collect();

    let mut matched_writer = matched_writer.lock().unwrap();
    let mut unmatched_writer = unmatched_writer.lock().unwrap();
    let mut matched_count = matched_count.lock().unwrap();
    let mut unmatched_count = unmatched_count.lock().unwrap();

    for (is_match, id, seq) in results {
        if is_match {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            *matched_count += 1;
        } else {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            *unmatched_count += 1;
        }
    }

    Ok(())
}

// Good print statements, second Rayon architecture
fn process_reads_rayon(
    reads_path: &str,
    processor: KmerProcessor,
    matched_path: &str,
    unmatched_path: &str,
) -> Result<(u32, u32, u32, u32), Box<dyn Error + Send + Sync>> {
    let mut reader = needletail::parse_fastx_file(reads_path)?;
    let processor = Arc::new(processor);

    let matched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(matched_path)?)));
    let unmatched_writer = Arc::new(Mutex::new(BufWriter::new(File::create(unmatched_path)?)));

    let mseq_count = Arc::new(Mutex::new(0u32));
    let mbase_count = Arc::new(Mutex::new(0u32));

    let useq_count = Arc::new(Mutex::new(0u32));
    let ubase_count = Arc::new(Mutex::new(0u32));

    let mut chunks = Vec::new();
    let mut chunk = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;
        let id = Arc::from(str::from_utf8(record.id())?);
        let seq = record.seq().into_owned();

        chunk.push((id, seq));

        if chunk.len() == 10_000 {
            chunks.push(chunk);
            chunk = Vec::new();
        }
    }

    if !chunk.is_empty() {
        chunks.push(chunk);
    }

    chunks.into_par_iter().try_for_each(|chunk| {
        process_chunk_rayon(
            chunk,
            processor.clone(),
            matched_writer.clone(),
            unmatched_writer.clone(),
            mseq_count.clone(),
            mbase_count.clone(),
            useq_count.clone(),
            ubase_count.clone(),
        )
    })?;

    let mseq_count = *mseq_count.lock().unwrap();
    let mbase_count = *mbase_count.lock().unwrap();

    let useq_count = *useq_count.lock().unwrap();
    let ubase_count = *ubase_count.lock().unwrap();

    Ok((mseq_count, mbase_count, useq_count, ubase_count))
}

fn process_chunk_rayon(
    reads_chunk: Vec<(Arc<str>, Vec<u8>)>,
    processor: Arc<KmerProcessor>,
    matched_writer: Arc<Mutex<BufWriter<File>>>,
    unmatched_writer: Arc<Mutex<BufWriter<File>>>,
    mseq_count: Arc<Mutex<u32>>,
    mbase_count: Arc<Mutex<u32>>,
    useq_count: Arc<Mutex<u32>>,
    ubase_count: Arc<Mutex<u32>>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let results: Vec<(bool, Arc<str>, Vec<u8>)> = reads_chunk
        .into_par_iter()
        .map(|(id, seq)| {
            let is_match = processor.process_read(&seq);
            (is_match, id, seq)
        })
        .collect();

    let mut matched_writer = matched_writer.lock().unwrap();
    let mut unmatched_writer = unmatched_writer.lock().unwrap();
    let mut mseq_count = mseq_count.lock().unwrap();
    let mut mbase_count = mbase_count.lock().unwrap();
    let mut useq_count = useq_count.lock().unwrap();
    let mut ubase_count = ubase_count.lock().unwrap();

    for (is_match, id, seq) in results {
        if is_match {
            writeln!(matched_writer, ">{}", id)?;
            matched_writer.write_all(&seq)?;
            writeln!(matched_writer)?;
            *mseq_count += 1;
            *mbase_count += seq.len() as u32;
        } else {
            writeln!(unmatched_writer, ">{}", id)?;
            unmatched_writer.write_all(&seq)?;
            writeln!(unmatched_writer)?;
            *useq_count += 1;
            *ubase_count += seq.len() as u32;
        }
    }

    Ok(())
}
