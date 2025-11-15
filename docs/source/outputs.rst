Output Files
============

Sheriff generates 13 output files in the specified output directory. This page describes each file in detail.

Output Directory Structure
---------------------------

.. code-block:: text

   output_directory/
   ├── edit_site_info.txt                                       # Edit site metadata
   ├── edit_sites.bed                                           # BED format
   ├── t7_barcode_edits.tsv                                     # Edit details
   ├── cell_allelic_dosage.canonical-edit-sites.parquet.gz      # Cell × edit counts
   ├── cell_allelic_dosage.canonical-edit-sites.gene-collapsed.parquet.gz
   ├── cell_gene_mrna_counts.parquet.gz                         # Gene expression
   ├── t7_barcoded_counts.parquet.gz                            # Barcoded UMI counts
   ├── t7_nonbarcoded_counts.parquet.gz                         # Non-barcoded UMI counts
   ├── t7_all_counts.parquet.gz                                 # All T7 UMIs
   ├── t7_reads.txt                                             # Read names
   ├── t7_barcoded_reads.txt
   ├── t7_non-barcoded_reads.txt
   ├── t7_filt.bam                                              # BAM subsets
   ├── t7_only.bam
   ├── t7_non-barcoded_only.bam
   └── t7_barcoded_only.bam

Primary Analysis Files
----------------------

edit_site_info.txt
^^^^^^^^^^^^^^^^^^

**Format:** Tab-separated text

**Description:** Comprehensive metadata for each canonical edit site.

**Columns:**

* ``edit_site``: Chromosome and position (e.g., "chr1:12345")
* ``chrom``: Chromosome name
* ``start``: Edit site start position
* ``end``: Edit site end position
* ``n_cells``: Number of cells with edits at this site
* ``forward_reads``: Count of forward-strand T7 reads
* ``reverse_reads``: Count of reverse-strand T7 reads
* ``stranded_edit_dist``: Distance between nearest forward/reverse edits
* ``overlapping_genes``: Comma-separated list of intersecting genes
* ``edit_type``: Classification (e.g., "canonical", "off-target")

**Usage:**

.. code-block:: python

   import pandas as pd
   edit_info = pd.read_csv("edit_site_info.txt", sep='\t')

   # Filter for high-confidence edits
   high_conf = edit_info[edit_info['n_cells'] > 10]

   # Get edits in specific gene
   gene_edits = edit_info[edit_info['overlapping_genes'].str.contains('BRCA1')]

edit_sites.bed
^^^^^^^^^^^^^^

**Format:** BED6

**Description:** Edit sites in standard BED format for genome browsers.

**Columns:**

1. Chromosome
2. Start position
3. End position
4. Name (edit_site identifier)
5. Score (number of cells)
6. Strand

**Usage:**

Load in IGV, UCSC Genome Browser, or other tools for visualization.

t7_barcode_edits.tsv
^^^^^^^^^^^^^^^^^^^^

**Format:** Tab-separated text

**Description:** Detailed information for individual edit events.

**Columns:**

* ``read_name``: BAM read name
* ``cell_barcode``: Cell barcode
* ``umi``: UMI sequence
* ``chrom``: Chromosome
* ``position``: Edit position
* ``ref_seq``: Reference sequence
* ``alt_seq``: Alternative sequence (with insertion)
* ``edit_site``: Associated canonical edit site
* ``kmer_matches``: Number of T7 k-mer matches

**Usage:**

.. code-block:: python

   edits = pd.read_csv("t7_barcode_edits.tsv", sep='\t')

   # Get unique edits per cell
   unique_edits = edits.groupby(['cell_barcode', 'edit_site']).size()

Count Matrices
--------------

cell_allelic_dosage.canonical-edit-sites.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Format:** Parquet (compressed)

**Dimensions:** cells × edit_sites

**Description:** UMI counts for edited alleles at each canonical edit site per cell.

**Values:** Integer counts (0 = no edit, >0 = number of edited UMI)

**Usage in Python:**

.. code-block:: python

   import pandas as pd
   allelic_dosage = pd.read_parquet(
       "cell_allelic_dosage.canonical-edit-sites.parquet.gz"
   )

   # Get cells with edit at specific site
   edit_site = "chr1:12345"
   cells_with_edit = allelic_dosage[allelic_dosage[edit_site] > 0]

   # Calculate edit frequencies
   edit_freq = (allelic_dosage > 0).sum(axis=0) / len(allelic_dosage)

**Usage in R:**

.. code-block:: r

   library(arrow)
   allelic_dosage <- as.data.frame(
       read_parquet("cell_allelic_dosage.canonical-edit-sites.parquet.gz")
   )
   rownames(allelic_dosage) <- allelic_dosage[, '__index_level_0__']
   allelic_dosage <- allelic_dosage[, 1:(ncol(allelic_dosage)-1)]

cell_allelic_dosage.canonical-edit-sites.gene-collapsed.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Format:** Parquet (compressed)

**Dimensions:** cells × (edit_sites + genes)

**Description:** Same as above, but edit sites overlapping genes are collapsed to the gene level.

**Use Case:** Simplifies analysis when multiple edits target the same gene.

cell_gene_mrna_counts.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Format:** Parquet (compressed)

**Dimensions:** cells × genes

**Description:** Gene expression UMI counts with T7 reads filtered out.

**Values:** Integer UMI counts

**Usage:**

.. code-block:: python

   gene_counts = pd.read_parquet("cell_gene_mrna_counts.parquet.gz")

   # Get expression of specific gene
   brca1_expr = gene_counts['BRCA1']

   # Cells expressing gene
   brca1_positive = gene_counts[gene_counts['BRCA1'] > 0]

   # Normalize to CPM
   cpm = gene_counts.div(gene_counts.sum(axis=1), axis=0) * 1e6

T7 Count Matrices
-----------------

t7_barcoded_counts.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** UMI counts for T7-barcoded reads at each edit site.

**Use:** High-confidence T7 reads with k-mer barcode matches.

t7_nonbarcoded_counts.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** UMI counts for reads near edit sites without T7 barcode.

**Use:** Potential T7 reads captured by proximity to canonical sites.

t7_all_counts.parquet.gz
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Sum of barcoded + non-barcoded T7 counts.

**Use:** Total T7-related UMIs per cell per edit site.

BAM Subsets
-----------

These files are useful for visualization and quality control:

t7_filt.bam
^^^^^^^^^^^

**Description:** BAM containing only mRNA reads (T7 reads removed).

**Use:** Standard gene expression analysis without T7 contamination.

t7_only.bam
^^^^^^^^^^^

**Description:** BAM containing all potential T7 reads.

**Use:** Inspect T7 read distribution in genome browser.

t7_barcoded_only.bam
^^^^^^^^^^^^^^^^^^^^^

**Description:** BAM containing only T7-barcoded reads.

**Use:** High-confidence T7 reads for visualization.

t7_non-barcoded_only.bam
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** BAM containing reads near edit sites without T7 barcode.

**Use:** Inspect non-barcoded T7 candidates.

Read Name Lists
---------------

t7_reads.txt
^^^^^^^^^^^^

**Description:** List of all T7 read names.

t7_barcoded_reads.txt
^^^^^^^^^^^^^^^^^^^^^

**Description:** List of T7-barcoded read names only.

t7_non-barcoded_reads.txt
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** List of non-barcoded T7 read names.

File Size Expectations
----------------------

For a typical 10k cell experiment:

.. list-table::
   :header-rows: 1

   * - File
     - Typical Size
   * - edit_site_info.txt
     - < 1 MB
   * - cell_allelic_dosage.*.parquet.gz
     - 10-100 MB
   * - cell_gene_mrna_counts.parquet.gz
     - 50-500 MB
   * - t7_filt.bam
     - Several GB
   * - t7_only.bam
     - 100-500 MB

Integration with Analysis Tools
--------------------------------

Scanpy (Python)
^^^^^^^^^^^^^^^

.. code-block:: python

   import scanpy as sc
   import pandas as pd

   # Load gene expression
   gene_counts = pd.read_parquet("cell_gene_mrna_counts.parquet.gz")
   adata = sc.AnnData(gene_counts)

   # Load edit dosage as obs
   edit_dosage = pd.read_parquet(
       "cell_allelic_dosage.canonical-edit-sites.parquet.gz"
   )
   adata.obs = pd.concat([adata.obs, edit_dosage], axis=1)

   # Standard analysis
   sc.pp.normalize_total(adata)
   sc.pp.log1p(adata)
   sc.pp.highly_variable_genes(adata)
   sc.tl.pca(adata)
   sc.tl.umap(adata)

Seurat (R)
^^^^^^^^^^

.. code-block:: r

   library(Seurat)
   library(arrow)

   # Load gene expression
   gene_counts <- as.data.frame(
       read_parquet("cell_gene_mrna_counts.parquet.gz")
   )
   rownames(gene_counts) <- gene_counts[, '__index_level_0__']
   gene_counts <- gene_counts[, -1]

   # Create Seurat object
   seurat_obj <- CreateSeuratObject(counts = t(gene_counts))

   # Load edit dosage
   edit_dosage <- as.data.frame(
       read_parquet("cell_allelic_dosage.canonical-edit-sites.parquet.gz")
   )
   rownames(edit_dosage) <- edit_dosage[, '__index_level_0__']
   edit_dosage <- edit_dosage[, -1]

   # Add to metadata
   seurat_obj <- AddMetaData(seurat_obj, metadata = edit_dosage)

See Also
--------

* :doc:`tutorial` - Example loading and analysis
* :doc:`usage` - Parameter descriptions
* :doc:`api/modules` - Programmatic access to outputs
