[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# **RustDUK**  
A high-performance Rust tool for filtering sequencing reads based on reference k-mers.
Inspired by [BBDuk](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) by Brian Bushnell. Provides performant and memory-efficient read processing with support for both paired and unpaired FASTA/FASTQ files, with multiple files or interleaved format.  

---

## **Features and Default Behavior**

- **K-mer based read filtering**:  
  - Reads are compared to reference sequences by matching k-mers.
  - If a read sequence has at least x k-mers also found in reference dataset, it is a match
  - x is 1 by default, changed with `--minhits <int>`

- **Paired reads support**:  
  - Paired inputs and outputs can be specified by adding more input/output files
  - Interleaved inputs or outputs, signify interleaved input with `--interinput`
  - Automatic detection of input/output mode

- **Multithreading with Rayon**:  
  - Adjustable thread count via `--threads` argument  
  - Defaults to all available cores

- **Memory Limit**:  
  - Specify maximum memory usage with `--maxmem <String>` (e.g., `5G` for 5 gigabytes, `500M` for 500 megabytes)  
  - Defaults to using 85% of available system memory

- **Automatic Reference Indexing**:  
  - Builds a serialized reference k-mer index using Bincode if `--binref <file>` is provided from references provided with `--ref <file>`
  - Uses saved index on subsequent runs if `--binref <file>` points to a serialized hashset of kmers
  
- **Output statistics**:
  - Total reads and bases processed  
  - Matches and non-matches  
  - Processing speed (reads/sec, bases/sec)

## **Performance Variables**
| **Parameter**        | **Default**                    | **Notes**                           |
|----------------------|--------------------------------|-------------------------------------|
| Threads              | Number of available threads    | Adjust with `--threads <int>`       |
| Max memory           | 85% of system memory           | Adjust with `--maxmem <String>`     |
| Chunk size           | 10,000 reads                   | Modify in source code               |
| Serialization format | [Bincode](https://docs.rs/bincode/) Standard | Modify in source code |

See more parameter documentation at ./rustduk.sh

---

## **Installation**

### **1. Install Rust**
If using UNIX, run this command and follow the ensuing instructions

`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`

If using Windows, download the correct installer from [Rustup](https://rustup.rs/#)

### **2. Clone the repository**
```bash
git clone https://github.com/jackdougle/rustduk.git
cd rustduk
```

### **3. Build the project**
Use the Rust release build for best performance:
```bash
cargo build --release
```

The binary will be located at:
```
target/release/rustduk
```

---

## **Example Usage**
```bash
./rustduk.sh --in reads.fq --ref reference.fa --outm matched.fq --outu unmatched.fq ...
```

This command:
- Reads input reads from `reads.fq`
- Filters them against reference k-mers built from `references.fa`
- Outputs matched reads to `matched.fq` and unmatched reads to `unmatched.fq`

---

## **Output**

Output files support either FASTA or FASTQ format.
- Will default to FASTQ unless extension is .fa, .fna, or .fasta

Example console output:

```
Indexing time:        0.040 seconds

Processing reads from reads.fq using 14 threads
Input and output is processed as interleaved
Processing time:      0.110 seconds

Input:                1000000 reads             150000000 bases
Matches:              20000 reads (2.00%)       3000000 bases (2.00%)
Nonmatches:           980000 reads (98.00%)     147000000 bases (98.00%)

Time:                 0.150 seconds
Reads Processed:      1.00m reads               6.65m reads/sec
Bases Processed:      150.00m bases             1000.00m bases/sec
```

---

### **Program Stages**

1. **Load reference k-mers**  
   - If serialized reference k-mer index exists, load using Bincode  
   - Else, parse reference FASTA and build k-mer index  
   - Serialize for future runs, if `--binref <File>` is provided

2. **Process reads in chunks**  
   - Reads are batched into chunks of 10,000 records  
   - Each chunk is then processed in parallel using Rayon

3. **Output matched/unmatched reads**  
   - Results are written as they are processed  
   - Separate or interleaved output depending on mode
   - Input second pair of output files for 2nd pair of paired reads output

4. **Statistics collection**  
   - Atomic counters track total reads and bases matched/unmatched  
   - Final summary printed at the end

### **Processing Modes**
RustDUK automatically detects the appropriate read handling mode:

| **Mode**                   | **Description**                                     |
|----------------------------|-----------------------------------------------------|
| `Unpaired`                 | Single input file, unpaired reads                   |
| `Paired`                   | Two input files, separate outputs                   |
| `PairedInInterOut`         | Two input files, interleaved output                 |
| `InterInPairedOut`         | Interleaved input file, separate paired outputs     |
| `Interleaved`              | Interleaved input and output                        |

---

## **Future Features**

- Piping via stdin/stdout for certain workflows  
- K-mer metadata for serialized k-mer file to prevent mismatched k-mer size

Request more by emailing jack.gdouglass@gmail.com

---

## **License**

This project is licensed under the MIT License, see the [LICENSE](LICENSE) file for details.  

---

## **Crates Used**

- [Needletail](https://github.com/onecodex/needletail) — FASTA/FASTQ file reading and parsing and bitkmer operations
- [Bincode](https://sr.ht/~stygianentity/bincode/) — K-mer hashset serialization/deserialization
- [Rayon](https://github.com/rayon-rs/rayon) — Multithreading
- [Rustc-Hash](https://github.com/rust-lang/rustc-hash) — FxHashset for k-mer storing and comparison
- [Clap](https://github.com/clap-rs/clap) — CLI
- [Num-Cpus](https://github.com/seanmonstar/num_cpus) — detection of available CPU cores
- [Sysinfo](https://github.com/GuillaumeGomez/sysinfo) — system memory and resource information
- [Bytesize](https://github.com/tailhook/bytesize) — human-readable byte size formatting