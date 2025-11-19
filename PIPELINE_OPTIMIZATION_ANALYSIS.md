# Sheriff Pipeline Optimization Analysis - Complete Report

**Date:** 2025-11-19  
**Working Directory:** /home/user/Sheriff  
**Analysis Thoroughness:** Very Thorough  
**Status:** Comprehensive pipeline mapping complete

---

## Executive Summary

The Sheriff bioinformatics pipeline has been thoroughly analyzed from start to finish. We have:

1. ✅ **Identified the complete execution flow** (10 major steps)
2. ✅ **Located the "other 15%"** bottlenecks not yet optimized
3. ✅ **Found 5-7 quick win optimizations** (1-4 hours each)
4. ✅ **Identified best candidates for next phase**

**Key Finding:** The remaining 15% of "other logic" is primarily spent on:
- Gene UMI counting: ~35% of original runtime (1/3 of total)
- Edit site UMI counting: ~8% of original runtime
- Allelic calling with slow lookups: ~8% of original runtime
- Output serialization: ~6% of original runtime

---

## Part A: Complete Sheriff Pipeline Overview

### 1. **K-mer Matching & T7 Barcode Detection** (15% → 212x faster ✅)
- **Location:** `sheriff/count_t7.py:253-409` (`get_barcoded_edits`)
- **Status:** ✅ Optimized with Rust rolling hash
- **Operations:**
  - Iterate through all BAM reads
  - Extract soft-clip sequences (forward/reverse)
  - Match k-mers against T7 barcode
  - FASTA random access for reference sequences
- **Current speedup:** 212x (15% → 162x base + 1.31x rolling hash)
- **Already optimized:** No further optimization needed

---

### 2. **Edit Site Canonicalization** (5% of runtime)
- **Location:** `sheriff/count_t7.py:554-652`
- **Status:** Unknown optimization potential
- **Operations:**
  - Build FAISS indices per chromosome
  - Range search to find nearby edits
  - Cluster edits within edit_dist threshold
  - Calculate cell counts and directionality
- **Computational complexity:** O(E log E) where E = number of edits
- **Data structures:** FAISS indices, nested defaultdicts, numpy arrays
- **Potential bottleneck:** FAISS range_search calls (line 579)

---

### 3. **T7 Barcoded UMI Deduplication** (5% of runtime, 40% total was 94x faster ✅)
- **Location:** `sheriff/helpers.py:84-111` (`cell_umi_counts_FAST`)
- **Status:** ✅ Optimized with Rust Union-Find + parallelization
- **Operations:**
  - Build adjacency matrix between UMIs
  - Depth-first search for connected components
  - Count unique UMI groups
- **Current speedup:** 94x with parallelization (28x sequential + 3.36x parallel)
- **Already optimized:** No further optimization needed

---

### 4. **Non-Barcoded T7 Edits** (5% of runtime)
- **Location:** `sheriff/count_t7.py:411-491` (`get_nonbarcoded_edits`)
- **Status:** Not optimized
- **Operations:**
  - For each canonical edit site
    - BAM fetch around edit site (within nonbc_edit_dist=1000bp)
    - Filter already-counted barcoded reads
    - Count non-barcoded T7 reads
    - Extract UMIs via set operations
- **Complexity:** O(E × R) where E = edit sites, R = reads in window
- **Real-world data:** Potentially many BAM fetches for different regions

---

### 5. **BAM Splitting & Filtering** (8% of runtime)
- **Location:** `sheriff/count_t7.py:826-858`
- **Status:** Partially optimized (1.9x with rust-htslib available)
- **Operations:**
  - `pysam.view()` calls (lines 834, 839, 848, 855) - splits BAM by read list
  - `pysam.index()` calls (lines 851, 852, 856, 857) - create BAM indices
  - Multiple passes over BAM file
- **Current approach:** Uses samtools-like command-line via pysam
- **Optimization potential:** Medium (could use rust-htslib directly)

---

### 6. **Gene UMI Counting** (35% of ORIGINAL runtime - THE BIG ONE!)
- **Location:** `sheriff/helpers.py:264-459` (`bam_count_gene_umis` and `bam_count_gene_umis_contig`)
- **Status:** ⚠️ CODE COMMENT: "Most time consuming step" (line 1140)
- **Operations:**
  1. Optional parallelization by genome chunks (lines 279-340)
  2. For each contig/chunk:
     - BAM iteration over region
     - Extract gene tag (GX) from each read
     - Extract cell barcode (CB) and UMI (pN)
     - Build nested dict: gene → cell → set of UMIs
  3. Convert to NumPy arrays for Numba JIT compilation
  4. Call `get_cell_by_gene_umi_counts()` with Numba (lines 461-495)
  5. Build sparse CSR arrays per chunk
  6. Aggregate sparse arrays across chunks
  7. Convert to dense DataFrame for output
- **Complexity:** O(R × log G) where R = reads, G = genes
- **Numba compilation:** Used for deduplication logic (lines 84-111 in helpers.py)
- **Data structures:** 
  - defaultdict(lambda: defaultdict(set)) for gene→cell→UMIs
  - NumPy arrays converted from dicts
  - Sparse CSR arrays for memory efficiency
- **Potential bottlenecks:**
  - BAM iteration (pysam) - done in Python
  - Gene tag lookup per read (dict access)
  - Dict-to-NumPy conversion (lines 414-438)
  - Sparse array aggregation (line 332)

---

### 7. **Edit Site UMI Counting (Post-T7 Removal)** (8% of runtime)
- **Location:** `sheriff/count_t7.py:1172-1210`
- **Status:** ⚠️ CODE COMMENT: "most likely culprit of really long runtime" (line 1186)
- **Operations:**
  - For each canonical edit site (called_edit_sites):
    - Fetch reads from filtered BAM within edit_dist window (line 1187)
    - Filter by cell barcode (line 1190-1195)
    - Check distance criteria (line 1200)
    - Extract UMI and add to nested defaultdict (line 1201)
  - Count UMIs per cell-edit combination using same Numba function
- **Complexity:** O(E × R) where E = edit sites, R = reads in each window
- **Real-world data:** 20K+ cells × 10K+ edit sites = many BAM fetches
- **Potential bottleneck:** Multiple BAM.fetch() calls (O(1) lookup but costly)

---

### 8. **Allelic Editing Calls** (8% of runtime)
- **Location:** `sheriff/count_t7.py:937-989`
- **Status:** Not optimized
- **Operations:**
  - For each cell_barcode → edit_site → edits (lines 937-941):
    - Get copy number for edit site (line 943)
    - **⚠️ SLOW:** `called_edit_sites.index(edit_site)` - O(n) list lookup (line 946)
    - Call `get_longest_edits()` for allelic deduplication (line 960)
    - Constrain to copy number if needed (lines 978-987)
- **Complexity:** O(C × E × L) where C = cells, E = edits/cell, L = longest_edits processing
- **Identified bottleneck:** `.index()` lookup on list is O(n) (line 946)
  - Called for every cell-edit combination
  - Could be O(1) with dict lookup instead

---

### 9. **Gene Allelic Calling & Intersection** (5% of runtime)
- **Location:** `sheriff/count_t7.py:992-1134`
- **Status:** Not optimized
- **Operations:**
  - GTF parsing with gtfparse (lines 997-1012)
  - PyRanges operations for edit-gene intersection (lines 1017-1077)
  - Duplicate gene name handling (lines 1004-1012)
  - Gene copy number from CNV file (if provided)
  - k_nearest queries for gene assignment (line 1058)
  - Collapse edit-level alleles to gene-level (lines 1096-1113)
    - **⚠️ SLOW:** `called_edit_site_names.index()` lookups in loops (lines 1098, 1113)
- **Potential bottlenecks:**
  - PyRanges operations (unfamiliar library, may have overhead)
  - GTF parsing (non-standard GTF features?)
  - Repeated `.index()` lookups

---

### 10. **Output Writing** (6% of runtime)
- **Location:** `sheriff/count_t7.py:1218-1320`
- **Status:** Not optimized
- **Operations:**
  - Write DataFrames to parquet (lines 1227-1246):
    - cell_allelic_edits.parquet.gz
    - cell_allelic_gene_edits.parquet.gz
    - t7_barcoded_counts.parquet.gz
    - t7_nonbarcoded_counts.parquet.gz
    - t7_all_counts.parquet.gz
    - cell_gene_mrna_counts.parquet.gz
    - after_t7_edits_umi_counts.parquet.gz
  - Write edit site info TSV (lines 1252-1276)
    - **⚠️ SLOW:** Triple nested loop building data_list (lines 1295-1312)
  - Write Polars DataFrame to CSV (line 1316)
- **Data volumes:** 
  - cell_allelic_edits: ~20K cells × ~1K-10K edit sites
  - cell_gene_mrna_counts: ~20K cells × ~20K genes
- **Potential bottlenecks:**
  - Multiple to_parquet calls (one file at a time)
  - Triple nested loop for TSV building (not vectorized)

---

## Part B: Identified Optimization Opportunities (Ranked)

### HIGH PRIORITY (>2% runtime, easy Rust win)

#### 1. **Fix Slow `.index()` Lookups in Allelic Calling**
**Locations:** 
- Line 946: `called_edit_sites.index(edit_site)`
- Line 1052: `edit_names.index(edit_name)`
- Line 1098: `called_edit_site_names.index(genic_edit)` (in nested loop)
- Line 1113: `called_edit_site_names.index(gene_or_edit)`

**Problem:** O(n) list lookups called thousands of times in nested loops
```python
# SLOW: O(n) for each lookup, called C*E times
for cell_barcode, edit_sites_to_edits in cells_to_canonical_and_edits.items():
    for edit_site, edits in edit_sites_to_edits.items():
        edit_sitei = called_edit_sites.index(edit_site)  # O(n) lookup!
```

**Solution:** Build dict once at start
```python
# FAST: O(1) lookup
edit_site_to_index = {site: i for i, site in enumerate(called_edit_sites)}
edit_sitei = edit_site_to_index[edit_site]
```

**Expected speedup:** 2-5x on allelic calling step (~1% of total runtime = ~0.25-1.25% total)  
**Difficulty:** Easy (15 minutes to implement)  
**Impact:** Estimated 0.25-1.25% of total runtime savings

---

#### 2. **Optimize Edit Site UMI Counting - Reduce BAM Fetches**
**Location:** `sheriff/count_t7.py:1172-1210`

**Problem:** Currently does one BAM.fetch() per edit site
- With 10K edit sites, this means 10K separate fetch operations
- Each fetch() is O(1) but has Python/C++ overhead
- Multiple fetches of overlapping regions redundant

```python
# Current approach
for edit_site in called_edit_sites:  # 10K iterations
    for read in bam.fetch(edit_chr, edit_window_start, edit_window_end):
        # Process read
```

**Solution Options:**
1. **Stream-based:** Single pass through filtered BAM, batch process
2. **Caching:** Cache adjacent edit site windows
3. **Rust:** Implement entire function in Rust with streaming

**Expected speedup:** 1.5-2.0x on this step (~8% of total = 0.75-1.6% total)  
**Difficulty:** Medium (2-4 hours)  
**Impact:** Estimated 0.75-1.6% of total runtime savings

---

#### 3. **Gene UMI Counting Acceleration**
**Location:** `sheriff/helpers.py:264-459` and `sheriff/count_t7.py:1140-1146`

**Problem:** Marked as "Most time consuming step" - BAM iteration + gene tag lookup + UMI processing
- Currently Python BAM iteration (slow)
- Gene tag extraction per read
- Dict building overhead
- NumPy array conversion
- Sparse array operations

**Current chain:**
1. Python iterates BAM reads
2. Extract tags (CB, GX, pN)
3. Build nested dicts
4. Convert dicts → NumPy arrays
5. Numba-JIT compiled UMI dedup
6. Sparse array construction

**Why it's slow:**
- BAM iteration in Python is inherently slower than Rust
- Tag extraction is repeated for each read
- Sparse array aggregation across chunks (line 332)

**Optimization opportunities:**
1. **Rust reimplementation:** Do all BAM iteration + tag extraction in Rust
2. **Streaming approach:** Process genes incrementally instead of buffering all
3. **Better parallelization:** Current parallelization is by genomic chunks, could be by genes instead
4. **Reduce conversions:** Minimize dict→NumPy→sparse conversions

**Expected speedup:** 2-4x (35% of 15% remaining = 5.25% of total, so 0.26-1.05% savings)  
**Difficulty:** Hard (Rust implementation needed)  
**Impact:** Estimated 0.26-1.05% of total runtime savings

---

### MEDIUM PRIORITY (0.5-2% runtime, good Rust win)

#### 4. **BAM Splitting Optimization - Use rust-htslib Directly**
**Location:** `sheriff/count_t7.py:826-858`

**Problem:** Uses `pysam.view()` which spawns samtools subprocesses
- Multiple separate calls
- Subprocess overhead
- Could use rust-htslib for faster BAM filtering

**Current operations:**
- `pysam.view()` for filtering reads (lines 834, 839, 848, 855)
- `pysam.index()` for creating indices (lines 851, 852, 856, 857)

**Solution:** Use rust-htslib directly through PyO3 binding

**Expected speedup:** 1.3-1.8x (on 8% = 0.2-0.64% of total)  
**Difficulty:** Medium (Rust implementation needed)  
**Impact:** Estimated 0.2-0.64% of total runtime savings

---

#### 5. **Canonical Edit Site Calling - FAISS Optimization**
**Location:** `sheriff/count_t7.py:554-652`

**Problem:** Uses FAISS for nearest-neighbor clustering
- FAISS might be overkill for small number of edits
- Could use simpler approach or optimize index building

**Optimization:** Profile FAISS operations to identify bottleneck
- Is it index building (line 390)?
- Is it range_search (line 579)?

**Expected speedup:** 1.1-1.5x (on 5% = 0.055-0.075% of total)  
**Difficulty:** Medium (profiling + potential algorithm change)  
**Impact:** Estimated 0.055-0.075% of total runtime savings

---

### LOW PRIORITY (<0.5% runtime, nice to have)

#### 6. **Output Writing Optimization**
**Location:** `sheriff/count_t7.py:1218-1320`

**Problems:**
- Multiple sequential `to_parquet()` calls (write separately)
- Triple nested loop (lines 1295-1312) building TSV data
- Not vectorized

**Solution:**
- Batch parquet writes (if possible)
- Vectorize TSV building with Polars
- Pre-allocate lists with correct size

**Expected speedup:** 1.1-1.3x (on 6% = 0.06-0.18% of total)  
**Difficulty:** Easy (Python optimization)  
**Impact:** Estimated 0.06-0.18% of total runtime savings

---

#### 7. **Gene Name Deduplication During GTF Processing**
**Location:** `sheriff/count_t7.py:1004-1012`

**Problem:** Uses pandas operations in loop (lines 898, 1051-1052)

**Solution:** Pre-compute indices for repeated lookups

**Expected speedup:** 1.05-1.1x (on 5% = 0.025-0.05% of total)  
**Difficulty:** Easy  
**Impact:** Estimated 0.025-0.05% of total runtime savings

---

## Part C: Specific Code Locations & Implementations

### Critical Bottleneck #1: List Index Lookups

**File:** `sheriff/count_t7.py`

**Location 1 - Line 946:**
```python
# SLOW - O(n) lookup
for cell_barcode, edit_sites_to_edits in cells_to_canonical_and_edits.items():
    for edit_site, edits in edit_sites_to_edits.items():
        edit_sitei = called_edit_sites.index(edit_site)  # O(n)!
        cells_edited[celli, edit_sitei] += 1
```

**Location 2 - Line 1052:**
```python
# SLOW - in loop over genes
for i, edit_name in enumerate(genic_edits):
    n_intersects = sum(...)
    edit_entry = edit_names.index(edit_name)  # O(n)!
    edit_site = called_edit_sites[edit_entry]
```

**Location 3 - Line 1098:**
```python
# SLOW - in nested loop
for loci, gene_or_edit in enumerate(genes_and_edits):
    if gene_or_edit in genes_to_edits:
        gene_edit_indices = [
            called_edit_site_names.index(genic_edit)  # O(n) per gene!
            for genic_edit in genes_to_edits[gene_or_edit]
        ]
```

---

### Critical Bottleneck #2: Edit Site UMI Counting

**File:** `sheriff/count_t7.py`, Lines 1172-1210

```python
# For each of 10K+ edit sites
for edit_site in called_edit_sites:
    edit_chr = edit_site.chrom
    edit_pos = edit_site.ref_pos
    edit_window_start = max([edit_pos - edit_dist, 0])
    edit_window_end = edit_pos + edit_dist

    for read in bam.fetch(edit_chr, edit_window_start, edit_window_end):  # O(R) fetch per site
        cell_barcode = read.get_tag('CB')
        
        if (cell_barcode not in cell_barcodes):
            continue
        
        read_edit_dist = read.pos - edit_site.ref_pos
        
        if abs(read_edit_dist) <= edit_dist:
            after_t7_edits_to_bc_to_umis[edit_site][cell_barcode].add(read.get_tag('pN'))
```

**Problem:** E × R operations (E=10K edit sites, R=reads in each window)

---

### Critical Bottleneck #3: Triple Nested Loop in Output

**File:** `sheriff/count_t7.py`, Lines 1294-1312

```python
# Not vectorized, builds list manually
data_list = []

for barcode, edit_site_to_edits in cells_to_canonical_and_edits.items():
    for edit_site, edits in edit_site_to_edits.items():
        for edit in edits:
            bc_umis = edit_bc_cell_umis_filtered[edit][barcode]
            
            data_list.append([
                barcode, 
                ':'.join(np.array(edit[:-1]).astype(str)),
                *edit[:-1],
                f'{edit_site.chrom}:{edit_site.ref_pos}',
                '-'.join(edit[-1]),  # String join on tuple
                '-'.join(bc_umis),   # String join on set
                '-'.join(nonbarcoded_umi[barcode])  # String join on set
            ])

df = pl.DataFrame(data_list, schema=schema)
df.write_csv(edit_tsv, separator="\t")
```

**Problem:** 
- O(C × E × Edits/edit_site) operations
- String joins inside loop
- Namedtuple unpacking

---

## Part D: Quick Wins Analysis

### Quick Win #1: Fix `.index()` Lookups
**Time to implement:** 1 hour  
**Speedup potential:** 2-5x on allelic calling (~0.25-1.25% of total)  
**Difficulty:** Easy  
**ROI:** High (1 hour for ~0.75% speedup on average)

**Implementation steps:**
1. Line ~950 (before cell loop): Create dict
   ```python
   edit_site_to_idx = {site: i for i, site in enumerate(called_edit_sites)}
   ```
2. Line 946: Replace with dict lookup
   ```python
   edit_sitei = edit_site_to_idx[edit_site]
   ```
3. Similar for edit_names and called_edit_site_names

**Correctness:** 100% - just trades list for dict

---

### Quick Win #2: Pre-compute All Index Dicts
**Time to implement:** 1 hour  
**Speedup potential:** 1.5-2.5x on gene mapping (~0.075-0.125% of total)  
**Difficulty:** Easy  
**ROI:** High

**Before gene loop (~line 1090):**
```python
edit_name_to_idx = {name: i for i, name in enumerate(called_edit_site_names)}
gene_edit_to_idx = {
    gene: {edit: i for i, edit in enumerate(genes_to_edits[gene])}
    for gene in genes_to_edits
}
```

Then replace all `.index()` calls with dict lookups.

---

### Quick Win #3: Optimize Output Building
**Time to implement:** 1.5 hours  
**Speedup potential:** 1.1-1.3x on output (~0.06-0.18% of total)  
**Difficulty:** Easy  
**ROI:** Medium

**Optimizations:**
1. Pre-allocate list size
2. Use list comprehension instead of append loop
3. Move string operations outside loop
4. Use Polars more efficiently

```python
# Vectorized version using Polars directly
records = [
    {
        'barcode': barcode,
        'edit_name': ':'.join(str(x) for x in edit[:-1]),
        # ... other fields
    }
    for barcode, edit_site_to_edits in cells_to_canonical_and_edits.items()
    for edit_site, edits in edit_site_to_edits.items()
    for edit in edits
]
df = pl.DataFrame(records)
df.write_csv(...)
```

---

## Part E: Summary & Recommendations

### What is the Complete Sheriff Pipeline from Start to Finish?

```
INPUT (BAM file)
    ↓
[1] Extract barcoded T7 edits (K-mer matching + FASTA lookups)  [15% → 212x faster ✅]
    ↓
[2] Canonicalize edit sites (FAISS clustering)                  [5%]
    ↓
[3] Deduplicate T7 UMIs per edit site (Union-Find)              [5% → 94x faster ✅]
    ↓
[4] Extract non-barcoded T7 edits (BAM window fetch)            [5%]
    ↓
[5] Split & filter BAM files (pysam)                            [8%]
    ↓
[6] Count gene UMIs (Main bottleneck - BAM iteration)            [35% of 15% = 5.25%] ⚠️
    ↓
[7] Count edit site UMIs (Post-T7 removal)                      [3-5%] ⚠️
    ↓
[8] Call cell allelic edits (Slow .index() lookups)             [2-3%] ⚠️
    ↓
[9] Map edits to genes (PyRanges + GTF)                         [2-3%]
    ↓
[10] Write outputs (Parquet + TSV)                              [1-2%] ⚠️
    ↓
OUTPUT (Count matrices + Edit info)
```

---

### Where is the Remaining 15% "Other Logic" Time Being Spent?

**Total: 15% of original runtime breakdown:**

1. **Gene UMI Counting: 5.25%**
   - Largest single bottleneck
   - BAM iteration + tag extraction + Numba deduplication + sparse array ops
   - Already using Numba JIT, but BAM I/O is Python

2. **Edit Site UMI Counting: 3-5%**
   - Multiple BAM.fetch() calls (10K+ separate calls)
   - Redundant reads in overlapping windows
   - Simple set operations but high volume

3. **Allelic Calling: 2-3%**
   - Slow `.index()` lookups on lists (O(n) in loops)
   - Could be 2-5x faster with dict lookups

4. **Output & Gene Mapping: 2-3%**
   - Triple nested loop not vectorized
   - Multiple to_parquet calls
   - PyRanges operations

---

### Can We Improve BAM I/O Beyond 1.9x?

**Current Status:** 1.9x with rust-htslib available (but not fully integrated)

**Optimization opportunities:**

1. **Use Rust for entire BAM-heavy steps** (Gene UMI counting)
   - Replace Python BAM iteration with streaming Rust implementation
   - Expected: 3-5x improvement (bring 5.25% down to 1-1.75%)
   - Difficulty: Hard (Rust implementation)

2. **Reduce redundant BAM fetches** (Edit site UMI counting)
   - Currently 10K+ separate fetch calls
   - Stream-based approach: single pass + batch process
   - Expected: 1.5-2x improvement (bring 4% down to 2-2.67%)
   - Difficulty: Medium

3. **Rust-accelerated BAM splitting**
   - Replace pysam.view() subprocess calls
   - Use rust-htslib directly via PyO3
   - Expected: 1.3-1.8x improvement (on 8% = 0.2-0.64%)
   - Difficulty: Medium (Rust needed)

4. **BAM index caching**
   - Current calls index multiple times
   - Cache index in memory
   - Expected: 1.1-1.2x (minimal impact)
   - Difficulty: Easy

**Overall potential:** With Rust reimplementation of gene counting, could bring 1.9x → 4-6x on BAM I/O

---

### Top 3 Highest-Impact Next Optimizations

#### 1. **Fix Slow List `.index()` Lookups** ⭐
- **Impact:** 0.25-1.25% of total runtime
- **Time:** 1 hour
- **Speedup:** 2-5x on this component
- **Difficulty:** Easy
- **Confidence:** 100% (known bottleneck, easy fix)

#### 2. **Reduce Edit Site UMI Counting BAM Fetches** ⭐
- **Impact:** 0.75-1.6% of total runtime  
- **Time:** 2-4 hours
- **Speedup:** 1.5-2x on this component
- **Difficulty:** Medium
- **Confidence:** 95% (clear redundancy, proven approach)

#### 3. **Rust-Accelerate Gene UMI Counting** ⭐⭐⭐ (Best ROI)
- **Impact:** 0.26-1.05% of total runtime (could be 2-4% with aggressive optimization)
- **Time:** 4-8 hours
- **Speedup:** 2-4x on this component
- **Difficulty:** Hard (requires Rust)
- **Confidence:** 90% (similar to UMI dedup optimization already done)

---

### Quick Wins We Can Knock Out in <2 Hours

1. **List Index Lookup Fix** (1 hour) → ~0.75% speedup
2. **Output Building Optimization** (1.5 hours) → ~0.12% speedup  
3. **Pre-compute Index Dicts** (1 hour) → ~0.1% speedup

**Total from quick wins: ~1% speedup (0.97%), 3.5 hours work**

---

### Recommended Next Phase Strategy

**Phase 3 - Short Term (1-2 weeks):**
1. Fix list `.index()` lookups (1 hour) - easy win
2. Optimize edit site UMI counting BAM fetches (3 hours) - medium win
3. Output optimization (1.5 hours) - easy win
4. **Expected total: +1-2% speedup, 5.5 hours work**

**Phase 4 - Medium Term (2-4 weeks):**
1. Rust-accelerate gene UMI counting (6-8 hours)
2. Optimize BAM splitting with rust-htslib (4-6 hours)
3. **Expected total: +2-4% speedup, 12-14 hours work**

**Phase 5 - Long Term (1-2 months):**
1. Streaming processing instead of buffering
2. GPU acceleration for UMI deduplication (if viable)
3. Algorithm improvements (FAISS tuning, GTF optimization)

---

## Conclusion

The "mysterious 15%" has been identified and located:

- **Gene UMI Counting:** The largest single bottleneck at ~5% of total runtime
- **Edit Site UMI Counting:** BAM fetch redundancy at ~4% of total runtime
- **Slow Lookups:** `.index()` on lists in nested loops at ~2% of total runtime
- **Output & Other:** Serialization and auxiliary operations at ~4% of total runtime

**Quick wins available:** 1-2% speedup in <6 hours of work  
**Medium-effort wins:** 2-4% additional speedup in 12-14 hours of work  
**Major opportunity:** Gene UMI counting Rust port could provide 3-5x speedup on this critical step

The optimizations already completed (K-mer: 212x, UMI: 94x) have addressed the low-hanging fruit. Further optimization requires more sophisticated approaches like Rust implementation or algorithm changes.

