Rust Acceleration API
=====================

Sheriff includes optional Rust acceleration that provides 10-100x speedup
for performance-critical operations while maintaining identical output.

Installation
------------

**Option 1: Install from source (recommended for development):**

.. code-block:: bash

   # Install Rust toolchain if not already installed
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env

   # Install maturin for building Rust-Python extensions
   pip install maturin

   # Build and install Rust acceleration
   cd sheriff-rs
   maturin develop --release
   cd ..

**Option 2: Future PyPI installation:**

.. code-block:: bash

   # Coming soon
   pip install sheriff[rust]

Overview
--------

Sheriff's Rust implementation provides dramatic performance improvements:

+------------------+-------------+-------------+----------+
| Operation        | Python Time | Rust Time   | Speedup  |
+==================+=============+=============+==========+
| BAM Filtering    | 180s        | 3.5s        | 51x      |
| K-mer Matching   | 60s         | 0.8s        | 75x      |
| Overall Pipeline | ~9 hours    | 1-2 hours   | 4.5-9x   |
+------------------+-------------+-------------+----------+

Automatic Fallback
~~~~~~~~~~~~~~~~~~

Sheriff automatically detects Rust availability:

- **Rust available:** Uses high-performance Rust implementations
- **Rust unavailable:** Falls back to pure Python (with performance warning)

This ensures Sheriff always works, regardless of Rust installation.

API Reference
-------------

BAM Filtering Functions
~~~~~~~~~~~~~~~~~~~~~~~

filter_bam_by_barcodes_rust
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   sheriff_rs.filter_bam_by_barcodes_rust(
       input_path: str,
       output_path: str,
       barcodes: List[str]
   ) -> Dict[str, Any]

Filter BAM file by cell barcode whitelist (sequential version).

**Arguments:**

- ``input_path`` (str): Path to input BAM file
- ``output_path`` (str): Path to output BAM file
- ``barcodes`` (List[str]): List of whitelisted cell barcodes

**Returns:**

Dictionary with:

- ``reads_processed`` (int): Total reads processed
- ``reads_kept`` (int): Reads matching whitelist
- ``reads_filtered`` (int): Reads removed
- ``duration_secs`` (float): Processing time in seconds

**Example:**

.. code-block:: python

   import sheriff_rs

   result = sheriff_rs.filter_bam_by_barcodes_rust(
       "input.bam",
       "filtered.bam",
       ["AAACCTGAGAAACCAT", "AAACCTGAGAAACCGC", ...]
   )

   print(f"Kept {result['reads_kept']:,} / {result['reads_processed']:,} reads")
   print(f"Duration: {result['duration_secs']:.2f}s")

filter_bam_by_barcodes_rust_parallel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   sheriff_rs.filter_bam_by_barcodes_rust_parallel(
       input_path: str,
       output_path: str,
       barcodes: List[str]
   ) -> Dict[str, Any]

Parallel BAM filtering using rayon (2-5x speedup over sequential).

**Arguments:** Same as ``filter_bam_by_barcodes_rust``

**Returns:** Same as ``filter_bam_by_barcodes_rust``

**Performance:**

- Small files (<1M reads): Similar to sequential (overhead dominates)
- Medium files (1-10M reads): 2-3x speedup
- Large files (>10M reads): 3-5x speedup

filter_bam_by_barcodes_rust_chromosome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   sheriff_rs.filter_bam_by_barcodes_rust_chromosome(
       input_path: str,
       output_path: str,
       barcodes: List[str],
       num_threads: Optional[int] = None
   ) -> Dict[str, Any]

Chromosome-based parallel BAM filtering (10-50x speedup on large files).

**Arguments:**

- ``input_path`` (str): Path to input BAM file (must be indexed .bam.bai)
- ``output_path`` (str): Path to output BAM file
- ``barcodes`` (List[str]): List of whitelisted cell barcodes
- ``num_threads`` (Optional[int]): Number of threads (default: auto-detect)

**Returns:** Same as ``filter_bam_by_barcodes_rust``

**Performance:**

- Small files (<1M reads): Slower due to overhead (~4s fixed cost)
- Large files (>100M reads): 10-20x speedup
- Very large files (>1B reads): 20-50x speedup

**Example:**

.. code-block:: python

   import sheriff_rs

   # Process 937M reads with chromosome parallelism
   result = sheriff_rs.filter_bam_by_barcodes_rust_chromosome(
       "large_dataset.bam",
       "filtered.bam",
       whitelist_barcodes,
       num_threads=16
   )

   # Expected: ~15-20x speedup vs Python (180s → 10s)

K-mer Matching Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

match_kmer_rust
^^^^^^^^^^^^^^^

.. code-block:: python

   sheriff_rs.match_kmer_rust(
       sequence: str,
       k: int,
       whitelist: Optional[List[int]] = None,
       output_hash: bool = True
   ) -> List[Union[int, str]]

Match k-mers in DNA sequence against whitelist.

**Arguments:**

- ``sequence`` (str): DNA sequence (A/C/G/T/N)
- ``k`` (int): K-mer length (typically 6-15 for barcode matching)
- ``whitelist`` (Optional[List[int]]): List of k-mer hashes to match against (None = all k-mers)
- ``output_hash`` (bool): Return hashes (True) or k-mer strings (False)

**Returns:**

List of matching k-mer hashes (int) or strings (str)

**Behavior:**

- Skips k-mers containing 'N' (ambiguous bases)
- Uses frequency array for efficient counting
- Filters by whitelist if provided
- Returns empty list if no matches

**Performance:**

- ~200 ns per sequence (vs 10+ μs in Python)
- 50-100x faster than recursive Python implementation
- Scales linearly with sequence length

**Example:**

.. code-block:: python

   import sheriff_rs

   # Find all k-mers in sequence
   matches = sheriff_rs.match_kmer_rust("AAACGTTT", k=4)
   # Returns: [1, 6, 27, 111, 191]

   # Match against T7 barcode whitelist
   whitelist = [27, 45, 109]  # Pre-computed k-mer hashes
   matches = sheriff_rs.match_kmer_rust(
       "AAACGTTT",
       k=4,
       whitelist=whitelist,
       output_hash=True
   )
   # Returns: [27]

   # Get k-mer strings instead of hashes
   matches_str = sheriff_rs.match_kmer_rust(
       "AAACGTTT",
       k=4,
       whitelist=[27],
       output_hash=False
   )
   # Returns: ["ACGT"]

**Typical Use in Sheriff:**

.. code-block:: python

   from sheriff.count_t7 import KmerMatcher, match_kmer

   # Create matcher with T7 barcode
   matcher = KmerMatcher(k=13, sequences=["GGGAGAGTAT"])

   # Match k-mers in indel sequence (automatically uses Rust if available)
   indel_seq = read.query_sequence[:read.query_alignment_start]
   matches = match_kmer(matcher, indel_seq, output_kmer_hash=True)

count_kmers_rust
^^^^^^^^^^^^^^^^

.. code-block:: python

   sheriff_rs.count_kmers_rust(
       sequence: str,
       k: int
   ) -> List[int]

Count all k-mer occurrences in sequence.

**Arguments:**

- ``sequence`` (str): DNA sequence
- ``k`` (int): K-mer length

**Returns:**

List of length 4^k with count for each possible k-mer

**Example:**

.. code-block:: python

   import sheriff_rs

   counts = sheriff_rs.count_kmers_rust("AAACGTTT", k=4)
   # Returns: [0, 1, 0, ..., 1, ...]  (length 256 for k=4)

Performance Tips
----------------

Choosing the Right BAM Filtering Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Sequential** (``filter_bam_by_barcodes_rust``):

   - Best for: Small files (<1M reads)
   - Minimal overhead
   - Simplest to use

2. **Parallel** (``filter_bam_by_barcodes_rust_parallel``):

   - Best for: Medium files (1-10M reads)
   - 2-5x speedup
   - No index required

3. **Chromosome** (``filter_bam_by_barcodes_rust_chromosome``):

   - Best for: Large files (>100M reads)
   - 10-50x speedup
   - Requires .bam.bai index
   - Most effective on production datasets

**Rule of thumb:**

- <1M reads → Sequential
- 1-10M reads → Parallel
- >10M reads → Chromosome (if indexed)

K-mer Matching Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Rust k-mer implementation is automatically used when available:

.. code-block:: python

   from sheriff.count_t7 import match_kmer, HAS_RUST_KMER

   # Check if Rust is available
   if HAS_RUST_KMER:
       print("✓ Using Rust k-mer matching (50-100x speedup)")
   else:
       print("✗ Using Python k-mer matching (slower)")

   # Function automatically chooses best implementation
   result = match_kmer(matcher, sequence, output_kmer_hash=True)

Memory Considerations
~~~~~~~~~~~~~~~~~~~~~

- BAM filtering: Memory usage ~1.5x vs Python (temporary chromosome files)
- K-mer matching: Memory usage identical to Python (4^k frequency array)

For very large datasets:

.. code-block:: python

   # Use chromosome mode with thread limit to control memory
   result = sheriff_rs.filter_bam_by_barcodes_rust_chromosome(
       input_bam,
       output_bam,
       barcodes,
       num_threads=8  # Limit concurrent chromosomes
   )

Building from Source
--------------------

For development or customization:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/BradBalderson/Sheriff.git
   cd Sheriff

   # Install Python dependencies
   pip install -e .

   # Build Rust module
   cd sheriff-rs
   cargo test  # Run Rust tests
   maturin develop --release  # Build and install
   cd ..

   # Run benchmarks
   python benchmarks/compare_all_modes.py

Troubleshooting
---------------

**"Rust k-mer matching not available" warning:**

Install Rust and build the module:

.. code-block:: bash

   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   source $HOME/.cargo/env
   cd sheriff-rs && maturin develop --release && cd ..

**"BAM file must be indexed" error:**

Create .bam.bai index:

.. code-block:: bash

   samtools index input.bam

**Slower than expected:**

- Small files have overhead (use sequential mode)
- Check if Rust module is installed correctly
- Ensure BAM file is indexed for chromosome mode

Citation
--------

If using Sheriff's Rust acceleration in publications, please cite:

  Lorenzini et al. (2025). Joint single-cell profiling of CRISPR-Cas9 edits
  and transcriptomes reveals widespread off-target events and their effects on
  gene expression. bioRxiv 2025.02.07.636966.
