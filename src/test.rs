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

fn parse_memory_size(size_str: &str) -> Result<u64, String> {
    let size_str = size_str.trim().to_uppercase();

    if let Some(num_str) = size_str.strip_suffix('G') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024 * 1024 * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else if let Some(num_str) = size_str.strip_suffix('M') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024 * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else if let Some(num_str) = size_str.strip_suffix('K') {
        num_str
            .parse::<u64>()
            .map(|n| n * 1024)
            .map_err(|_| format!("Invalid number: {}", num_str))
    } else {
        size_str
            .parse::<u64>()
            .map_err(|_| format!("Invalid size format: {}", size_str))
    }
}

#[cfg(unix)]
fn try_set_memory_limit(max_bytes: u64) -> io::Result<()> {
    // First try RLIMIT_RSS (physical memory) - more reliable on macOS
    if try_set_limit(RLIMIT_RSS, max_bytes, "physical").is_ok() {
        return Ok(());
    }

    // If RSS fails, try RLIMIT_AS (virtual memory) - works better on Linux
    if try_set_limit(RLIMIT_AS, max_bytes, "virtual").is_ok() {
        return Ok(());
    }

    Err(io::Error::new(
        io::ErrorKind::Unsupported,
        "System does not support rlimit memory limiting",
    ))
}

#[cfg(unix)]
fn try_set_limit(resource: i32, max_bytes: u64, mem_type: &str) -> io::Result<()> {
    // Get current limits
    let mut current_limit = rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };

    let res = unsafe { getrlimit(resource, &mut current_limit) };
    if res != 0 {
        return Err(io::Error::last_os_error());
    }

    if current_limit.rlim_cur != RLIM_INFINITY && current_limit.rlim_cur <= max_bytes {
        println!(
            "Memory limit already set to {:.1} MB via {}",
            current_limit.rlim_cur as f64 / 1024.0 / 1024.0,
            mem_type
        );
        return Ok(());
    }

    // Set the new limit
    let new_limit = rlimit {
        rlim_cur: max_bytes,
        rlim_max: current_limit.rlim_max,
    };

    println!(
        "Attempting to set {} memory limit to {:.1} MB",
        mem_type,
        max_bytes as f64 / 1024.0 / 1024.0
    );

    let res = unsafe { setrlimit(resource, &new_limit) };
    if res != 0 {
        let err = io::Error::last_os_error();
        return Err(err);
    }

    println!(
        "Successfully set {} limit to {:.1} MB",
        mem_type,
        max_bytes as f64 / 1024.0 / 1024.0
    );
    Ok(())
}

#[cfg(not(unix))]
fn try_set_memory_limit(_max_bytes: u64) -> io::Result<()> {
    println!("Memory limiting only works on UNIX systems\nRunning with default memory allocation");
    Ok(()) // Don't fail, just warn
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    if let Some(ref mem_str) = args.maxmem {
        match parse_memory_size(mem_str) {
            Ok(max_memory) => match try_set_memory_limit(max_memory) {
                Ok(()) => {}
                Err(e) => {
                    eprintln!("Could not set memory limit: {}", e);
                    eprintln!("Running with default memory allocation");
                }
            },
            Err(e) => {
                eprintln!("Error parsing memory size '{}': {}", mem_str, e);
                eprintln!(
                    "Valid formats: 'XG', 'XM', 'XK', or number with no suffix representing bytes"
                );
                std::process::exit(1);
            }
        }
    }

    duk::run(args)?;

    Ok(())
}
