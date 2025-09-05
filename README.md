[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# **Rust-DUK**  
A high-performance Rust tool for filtering sequencing reads based on reference k-mers, inspired by [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) by Brian Bushnell.  
Rust-DUK provides fast, memory-efficient read processing with support for both paired and unpaired FASTA/FASTQ files.  

---

## **Features**

- **K-mer based read filtering**:  
  - Reads are compared to reference sequences by matching k-mers.
  - If a read sequence has at least x hits, it is a match
  - x is 1 by default, changed with `--minhits <int>`
- **Paired reads support**:  
  - Paired-end inputs and outputs  
  - Interleaved inputs or outputs, signify interleaved input with `--interinput`
  - Automatic detection of input/output mode
- **Multithreading with Rayon**:  
  - Adjustable thread count via `--threads` argument  
  - Defaults to all available cores
- **Memory limit control**:  
  - Specify maximum memory usage with `--maxmem <String>` (e.g., `5G`, `500M`)  
  - Defaults to **85% of system memory**
- **Automatic reference indexing**:  
  - Builds a serialized reference k-mer index on first run  
  - Automatically reloads saved index on subsequent runs
- **Streaming support** (future):  
  Planned support for piping input/output via stdin/stdout.
- **Detailed statistics**:
  - Total reads and bases processed  
  - Matches and non-matches  
  - Processing speed (reads/sec, bases/sec)

---

## **Installation**

### **1. Install Rust**
If using UNIX, run this command and follow the ensuing instructions
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
If using Windows, download the correct installer from Rustup
`https://rustup.rs/#`

### **2. Clone the repository**
```bash
git clone https://github.com/jackdougle/rust-duk.git
cd rust-duk
```

### **3 Build the project**
Use the Rust release build for best performance:
```bash
cargo build --release
```

The binary will be located at:
```
target/release/rust-duk
```

---

## **Usage**

### **Basic command**
```bash
rust-duk --in reads.fq --ref reference.fna --outm matched.fq --outu unmatched.fq
```

This command:
- Reads input reads from `reads.fq`
- Filters them against reference k-mers built from `references.fa`
- Outputs matched reads to `matched.fq` and unmatched reads to `unmatched.fq`

---

### **Paired-end reads**
Provide two input files for paired-end data:
```bash
rust-duk \
  --in reads_R1.fq \
  --in2 reads_R2.fq \
  --outm matched_R1.fq \
  --outm2 matched_R2.fq \
  --outu unmatched_R1.fq \
  --outu2 unmatched_R2.fq
```

---

### **Interleaved input**
If your paired-end reads are interleaved into a single file:
```bash
rust-duk \
  --in reads_interleaved.fq \
  --outm matched.fq \
  --outu unmatched.fq \
  --interleaved
```

---

### **Thread control**
By default, Rust-DUK uses all available logical cores.  
You can limit the number of threads with `--threads`:
```bash
rust-duk \
  --in reads.fq \
  --ref reference.fna \
  --outm matched.fq \
  --outu unmatched.fq \
  --threads 8
```

---

### **Memory limit**
By default, Rust-DUK uses **85% of available system memory**.  
To manually specify a memory cap:
```bash
rust-duk \
  --in reads.fq \
  --ref reference.fna \
  --outm matched.fq \
  --outu unmatched.fq \
  --maxmem 5G
```

---

### **Serialized reference index**
The first time Rust-DUK is run with a new reference, it builds a binary k-mer index and saves it for faster loading in future runs.

Default behavior:
1. Looks for a serialized index at the path specified by `--binref`.
2. If not found, it will:
   - Load k-mers directly from the reference FASTA/FASTQ (`--ref`)
   - Serialize and save them for future runs to path specified with `--binref`

Example:
```bash
rust-duk \
  --in reads.fq \
  --ref references.fa \
  --binref reference_serialized.bin \
  --outm matched.fq \
  --outu unmatched.fq
```

---

## **Output**

Rust-DUK reports statistics after processing, for example:

```
Indexing time:		0.041 seconds

Processing reads from reads.fq using 14 threads
Input and output is processed as interleaved
Processing time:	0.109 seconds

Input:			        1000000 reads			    150000000 bases
Matches:		        20000 reads (2.00%) 		3000000 bases (2.00%)
Nonmatches:		        980000 reads (98.00%)		147000000 bases (98.00%)

Time:			        0.150 seconds
Reads Processed:	    1.00m reads			        6.65m reads/sec
Bases Processed:	    150.00m bases			    1000.00m bases/sec
```

---

## **Performance Tuning**

| **Parameter**        | **Default**                    | **Notes**                           |
|----------------------|--------------------------------|-------------------------------------|
| Threads              | Number of logical CPU cores    | Adjust with `--threads <int>`       |
| Max memory           | 85% of system memory           | Adjust with `--maxmem <String>`     |
| Chunk size           | 10,000 reads                   | Modify in source code               |
| Serialization format | [Bincode](https://docs.rs/bincode/) | Modify in source code               |

---

## **Anatomy**

### **Processing Modes**
Rust-DUK automatically detects the appropriate read handling mode:

| **Mode**                   | **Description**                                     |
|----------------------------|-----------------------------------------------------|
| `Unpaired`                 | Single input file, unpaired reads                   |
| `Paired`                   | Two input files, separate outputs                   |
| `PairedInInterOut`         | Two input files, interleaved output                 |
| `InterInPairedOut`         | Interleaved input file, separate paired outputs     |
| `Interleaved`              | Interleaved input and output                        |

---

### **Program Stages**

1. **Load reference k-mers**  
   - If serialized index exists: load via Bincode  
   - Else, parse reference FASTA and build index  
   - Serialize for future runs

2. **Process reads in chunks**  
   - Reads are batched into chunks of 10,000 records  
   - Each chunk is processed in parallel using Rayon

3. **Output matched/unmatched reads**  
   - Results are written as they are processed  
   - Separate or interleaved output depending on mode
   - Input second pair of output files for alternate paired read output

4. **Statistics collection**  
   - Atomic counters track total reads and bases matched/unmatched  
   - Final summary printed at the end

---

### **Memory Management**
- By default, Rust-DUK dynamically allocates memory up to **85% of total system RAM**.
- If `--maxmem` is provided, memory usage is capped accordingly.
- Memory is used primarily for:
  - Reference k-mer index
  - Input/output buffers
  - Parallel processing of read chunks

---

### **Thread Pool Management**
Rust-DUK uses Rayon for multithreaded processing
- Number of threads used can be controlled via `--threads`
- Defaults to all avilable threads

---

## **Future Features**

- [ ] Streaming via stdin/stdout for piping workflows  
- [ ] Adaptive chunk sizing for large datasets

Request more by emailing jack.gdouglass@gmail.com

---

## **License**

This project is licensed under the MIT License, see the [LICENSE](LICENSE) file for details.  

---

## **Crates Used**

- [Needletail](https://github.com/onecodex/needletail) — FASTA/FASTQ file reading and parsing and bitkmer operations
- [Bincode](https://sr.ht/~stygianentity/bincode/) — efficient serialization/deserialization
- [Rayon](https://github.com/rayon-rs/rayon) — multithreaded processing
- [Rustc-Hash](https://github.com/rust-lang/rustc-hash) — fast hash maps and sets for k-mer indexing
- [Clap](https://github.com/clap-rs/clap) — command line parsing
- [Num-Cpus](https://github.com/seanmonstar/num_cpus) — detection of available CPU cores
- [Sysinfo](https://github.com/GuillaumeGomez/sysinfo) — system memory and resource information
- [Bytesize](https://github.com/tailhook/bytesize) — human-readable byte size formatting

- Inspired by [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) from the BBTools suite of bioinformatics tools
