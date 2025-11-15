Usage Guide
===========

This guide covers all Sheriff command-line parameters and common usage patterns.

Basic Command Structure
-----------------------

.. code-block:: bash

   sheriff <bam_file> <ref_file> <barcode_file> <gtf_file> [OPTIONS]

Required Arguments
------------------

BAM File
^^^^^^^^

**Type:** Path to BAM file

The aligned and annotated BAM file from split-pipe (Parse Biosciences pipeline).

**Requirements:**

* Must be sorted and indexed (.bai file present)
* Must contain cell barcode (CB) tags
* Must contain UMI (UB) tags
* Output from split-pipe pipeline

Reference Genome
^^^^^^^^^^^^^^^^

**Type:** Path to FASTA file

The reference genome in FASTA format.

**Requirements:**

* Must be indexed (.fai file present)
* Should match the genome used for alignment
* Typically GRCh38 for human data

.. code-block:: bash

   samtools faidx reference.fa

Barcode Whitelist
^^^^^^^^^^^^^^^^^

**Type:** Path to text file

File containing whitelisted cell barcodes, one per line.

**Format:**

.. code-block:: text

   AAACCTGAGAAACCAT
   AAACCTGAGAAACCGC
   AAACCTGAGAAACCTG
   ...

GTF Annotation
^^^^^^^^^^^^^^

**Type:** Path to GTF file

Gene annotation file in GTF format.

**Source:** Typically from Ensembl or GENCODE

Optional Parameters
-------------------

T7 Barcode Detection
^^^^^^^^^^^^^^^^^^^^

``--t7 <sequence>``

Target barcode sequence to identify T7 reads.

**Default:** ``GGGAGAGTAT``

.. code-block:: bash

   sheriff ... --t7 GGGAGAGTAT

``-k, --kmer <size>``

Size of k-mers for barcode matching.

**Default:** 6

.. code-block:: bash

   sheriff ... -k 6

Edit Site Calling
^^^^^^^^^^^^^^^^^

``--edit_dist <distance>``

Distance threshold for grouping edits as the same site.

**Default:** 140

.. code-block:: bash

   sheriff ... --edit_dist 140

``--edit_site_min_cells <count>``

Minimum cells required to call a canonical edit site.

**Default:** 3

**Recommendation:**

* Use 3 for 10k cell libraries
* Use 1 for smaller libraries or exploration

.. code-block:: bash

   sheriff ... --edit_site_min_cells 3

``--bidirectional_inserts / --no-bidirectional_inserts``

Require bi-directional donor insertion evidence.

**Default:** Enabled

**Highly recommended** to keep enabled to reduce false positives.

.. code-block:: bash

   # Enable (default)
   sheriff ... --bidirectional_inserts

   # Disable (not recommended)
   sheriff ... --no-bidirectional_inserts

``--stranded_edit_dist <distance>``

Maximum allowed distance between nearest forward and reverse edit sites.

**Default:** 15

Only applies when ``--bidirectional_inserts`` is enabled.

.. code-block:: bash

   sheriff ... --stranded_edit_dist 15

``--nonbc_edit_dist <distance>``

Distance around edit sites to capture non-barcoded reads.

**Default:** 1000

.. code-block:: bash

   sheriff ... --nonbc_edit_dist 1000

Filtering Options
^^^^^^^^^^^^^^^^^

``--blacklist <bed_file>``

BED file specifying genomic regions to exclude.

**Use:** Filter out regions with high endogenous T7-like sequences that cause false positives.

.. code-block:: bash

   sheriff ... --blacklist blacklist_regions.bed

``--whitelist <bed_file>``

BED file specifying known edit sites.

**Use:** Regions intersecting these will be called as canonical edit sites regardless of other filters.

.. code-block:: bash

   sheriff ... --whitelist known_edit_sites.bed

``--blacklist_seqs <text_file>``

Text file of sequences to filter from soft-clip analysis (one per line).

**Use:** Filter common artifacts like TSO sequences.

.. code-block:: bash

   sheriff ... --blacklist_seqs artifact_seqs.txt

Ploidy and Copy Number
^^^^^^^^^^^^^^^^^^^^^^

``--ploidy <number>``

Expected number of chromosome copies.

**Default:** 2 (diploid)

.. code-block:: bash

   sheriff ... --ploidy 2

``--cnv <bedgraph_file>``

BedGraph file specifying copy number variants.

**Use:** Adjust expected allele counts at CNV regions.

.. code-block:: bash

   sheriff ... --cnv K562_CNV.tsv

Gene Counting Options
^^^^^^^^^^^^^^^^^^^^^

``--mrna_count_mode <mode>``

Mode for quantifying gene expression.

**Options:**

* ``all``: Count all reads (default)
* ``polyT``: Count only polyT reads (mature mRNA)

.. code-block:: bash

   sheriff ... --mrna_count_mode all

``--uncorrected_count / --no-uncorrected_count``

Output gene counts WITHOUT filtering T7 reads.

**Default:** Disabled

**Use:** For comparison with standard scRNA-seq pipelines.

.. code-block:: bash

   sheriff ... --uncorrected_count

Performance Options
^^^^^^^^^^^^^^^^^^^

``--cpu <number>``

Number of CPUs for parallel processing.

**Default:** 1

**Recommendation:** Use 4-8 CPUs for faster processing.

.. code-block:: bash

   sheriff ... --cpu 8

``--chunk <megabases>``

Genome chunk size in megabases for gene counting.

**Default:** 250

**Use:** Lower if encountering memory issues.

.. code-block:: bash

   sheriff ... --chunk 250

Output Options
^^^^^^^^^^^^^^

``-o, --outdir <directory>``

Output directory for results.

**Default:** Current working directory

.. code-block:: bash

   sheriff ... -o ./results/

``-v, --verbosity <level>``

Logging verbosity level.

**Options:**

* ``0``: Errors only
* ``1``: Progress messages (default)
* ``2``: Debug information

.. code-block:: bash

   sheriff ... -v 2

Common Usage Patterns
---------------------

Standard Analysis
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sheriff aligned.bam reference.fa barcodes.txt genes.gtf \
       --cpu 8 \
       --edit_site_min_cells 3 \
       --blacklist blacklist.bed \
       -o results/

With Copy Number Variants
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sheriff aligned.bam reference.fa barcodes.txt genes.gtf \
       --cpu 8 \
       --cnv cnv_regions.tsv \
       --ploidy 2 \
       -o results/

Exploration Mode (Permissive)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sheriff aligned.bam reference.fa barcodes.txt genes.gtf \
       --edit_site_min_cells 1 \
       --no-bidirectional_inserts \
       -v 2 \
       -o exploration_results/

High-Stringency Mode
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sheriff aligned.bam reference.fa barcodes.txt genes.gtf \
       --edit_site_min_cells 5 \
       --stranded_edit_dist 10 \
       --blacklist blacklist.bed \
       --whitelist known_sites.bed \
       -o high_stringency_results/

Tips and Best Practices
------------------------

**Performance:**

* Use ``--cpu`` to parallelize gene UMI counting
* Increase ``--chunk`` on high-memory systems
* Process smaller cell subsets for testing

**Quality Control:**

* Always use blacklist files to filter problematic regions
* Use ``--bidirectional_inserts`` (default) to reduce false positives
* Set ``--edit_site_min_cells`` based on library size (3 for 10k cells)

**Troubleshooting:**

* Increase ``-v`` to 2 for debugging
* Check BAM file has CB and UB tags
* Verify reference genome matches alignment
* Ensure all input files are properly indexed

See Also
--------

* :doc:`tutorial` - Step-by-step example
* :doc:`outputs` - Output file descriptions
* :doc:`troubleshooting` - Common issues
