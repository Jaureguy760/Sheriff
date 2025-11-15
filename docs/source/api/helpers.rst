helpers Module
==============

Helper functions for UMI counting and data processing.

.. automodule:: sheriff.helpers
   :members:
   :undoc-members:
   :show-inheritance:

UMI Counting Functions
----------------------

get_t7_count_matrix
^^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.get_t7_count_matrix

Builds cell × edit-site count matrices from UMI data.

get_cell_counts_from_umi_dict
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.get_cell_counts_from_umi_dict

Aggregates UMI counts per cell with deduplication.

bam_count_gene_umis
^^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.bam_count_gene_umis

Parallelized gene expression quantification.

**Features:**

* Multi-threaded processing
* Sparse matrix output for memory efficiency
* Chromosome-based chunking

UMI Deduplication
-----------------

deduplicate_umis
^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.deduplicated_umi_count_FAST

Collapses UMIs with single mismatches using graph-based approach.

**Algorithm:**

1. Build adjacency matrix of UMIs within 1 edit distance
2. Identify connected components (PCR duplicates)
3. Count unique UMI families

**Performance:** Optimized with Numba JIT compilation.

Filtering Functions
-------------------

get_longest_edits
^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.get_longest_edits

Identifies most significant edits per site based on insertion length.

bed_file_flag_edits
^^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.helpers.bed_file_flag_edits

Applies whitelist/blacklist filters from BED files.

Performance Notes
-----------------

Many functions in this module are optimized with Numba's ``@jit(nopython=True)`` decorator for high performance:

* ``cell_umi_counts_FAST``: Fast UMI counting
* ``deduplicated_umi_count_FAST``: Fast deduplication
* ``get_adjacency_matrix``: Adjacency matrix construction
* ``within_single_mismatch_int``: Edit distance calculation
* ``depth_first_search``: Graph traversal

**Memory Efficiency:**

* Uses ``scipy.sparse.csr_array`` for large count matrices
* Uses ``np.uint32`` for count data to reduce memory footprint
* Processes data in chunks to avoid loading entire genome

Examples
--------

**Basic UMI counting:**

.. code-block:: python

   from sheriff.helpers import get_cell_counts_from_umi_dict
   import numpy as np

   # Cell barcode to UMI mapping
   cell_bc_to_umis = {
       'AAACCTGAGAAACCAT': {'ATCGATCG', 'ATCGATCG', 'GCTAGCTA'},
       'AAACCTGAGAAACCGC': {'TTAACCGG', 'TTAACCGG'},
   }

   # Cell barcode index mapping
   cell_barcodes_dict = {
       'AAACCTGAGAAACCAT': 0,
       'AAACCTGAGAAACCGC': 1,
   }

   # Count unique UMIs per cell
   counts = get_cell_counts_from_umi_dict(cell_bc_to_umis, cell_barcodes_dict)
   # Returns: array([2, 1], dtype=uint32)

**Gene UMI counting from BAM:**

.. code-block:: python

   from sheriff.helpers import bam_count_gene_umis
   import pysam

   bam = pysam.AlignmentFile("data.bam")
   cell_barcodes = ['AAACCTGAGAAACCAT', 'AAACCTGAGAAACCGC']

   # Count gene UMIs (parallelized)
   gene_counts = bam_count_gene_umis(
       bam=bam,
       gtf_file="genes.gtf",
       cell_barcodes_list=cell_barcodes,
       n_cpus=4,
       chunk_size_mb=250,
       verbosity=1
   )
