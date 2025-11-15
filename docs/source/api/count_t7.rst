count_t7 Module
===============

Core processing module for T7 edit detection and counting.

.. automodule:: sheriff.count_t7
   :members:
   :undoc-members:
   :show-inheritance:

Main Function
-------------

run_count_t7
^^^^^^^^^^^^

.. autofunction:: sheriff.count_t7.run_count_t7

This is the main entry point for Sheriff processing. All command-line arguments are exposed as function parameters.

**Example:**

.. code-block:: python

   from sheriff.count_t7 import run_count_t7

   run_count_t7(
       bam_file="data.bam",
       ref_file="genome.fa",
       barcode_file="cells.txt",
       gtf_file="genes.gtf",
       t7_barcode="GGGAGAGTAT",
       k=6,
       edit_dist=140,
       stranded_edit_dist=15,
       edit_site_min_cells=3,
       nonbc_edit_dist=1000,
       ploidy=2,
       copy_number_variant_file=None,
       blacklist_file=None,
       whitelist_file=None,
       blacklist_seqs=None,
       mrna_count_mode="all",
       outdir="./output/",
       edit_site_rev_comp_filt=True,
       uncorrected_gene_count=False,
       verbosity=1,
       n_cpus=4,
       chunk_size_mb=250
   )

Helper Functions
----------------

get_barcoded_edits
^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.count_t7.get_barcoded_edits

Identifies T7-barcoded reads indicating edit sites.

get_nonbarcoded_edits
^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: sheriff.count_t7.get_nonbarcoded_edits

Captures non-barcoded reads near edit sites.

Classes
-------

KmerMatcher
^^^^^^^^^^^

.. autoclass:: sheriff.count_t7.KmerMatcher
   :members:
   :undoc-members:

Implements k-mer hashing and matching for barcode detection.

**Example:**

.. code-block:: python

   from sheriff.count_t7 import KmerMatcher

   # Create k-mer matcher
   matcher = KmerMatcher(k=6, sequences=["GGGAGAGTAT"])

   # Match k-mers in a sequence
   seq = "ATCGGGGAGAGTATATCG"
   matches = matcher.match_kmer(seq)
