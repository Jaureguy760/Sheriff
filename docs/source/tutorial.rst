Tutorial
========

This tutorial walks you through a complete Sheriff analysis using the included example data.

Overview
--------

**Expected Runtime:** < 30 seconds

**What You'll Learn:**

* How to prepare reference genome and annotations
* How to run Sheriff with appropriate parameters
* How to interpret the output files
* How to load results in Python or R

Prerequisites
-------------

* Sheriff installed (see :doc:`installation`)
* ~2 minutes for data preparation
* Basic familiarity with genomic data formats

Step 1: Prepare Reference Genome
---------------------------------

Download and index the human reference genome:

.. code-block:: bash

   # Navigate to Sheriff directory
   cd Sheriff

   # Download human genome (GRCh38)
   wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
       -O example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

   # Decompress
   gzip -d example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

   # Index with samtools
   mamba install samtools
   samtools faidx example_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa

Step 2: Prepare Gene Annotations
---------------------------------

Download the GTF annotation file:

.. code-block:: bash

   # Download GTF
   wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz \
       -O example_data/Homo_sapiens.GRCh38.110.gtf.gz

   # Decompress
   gzip -d example_data/Homo_sapiens.GRCh38.110.gtf.gz

Step 3: Set Up Parameters
--------------------------

The example data is from a 500-cell library with reads filtered to 200kb regions around true edit sites:

.. code-block:: bash

   # Define input files
   dir_="example_data/"
   bam_="${dir_}barcode_headAligned_anno.sorted.edit_regions_200kb.bam"
   ref_="${dir_}Homo_sapiens.GRCh38.dna.primary_assembly.fa"
   cells_="${dir_}barcode_whitelist.500-cell.txt"
   gtf_="${dir_}Homo_sapiens.GRCh38.110.gtf"

   # Optional filtering files
   cnv_="${dir_}K562_CNV_hg19.tsv"
   blacklist_="${dir_}black_100x_peaks_by_qval.simple_repeats_50N.EXTRA.bed"
   blacklist_seqs="${dir_}blacklist_seqs.txt"

   # Parameters
   min_="1"  # Minimum cells per edit site (use 3 for full dataset)
   cpu="1"   # Number of CPUs
   out_dir="./subset_500_cell_sheriff_output/"

.. note::
   For the 10k cell library, use ``--edit_site_min_cells 3``. The CNV file is optional and only needed if your cell line has copy number variants.

Step 4: Run Sheriff
--------------------

Execute Sheriff with the configured parameters:

.. code-block:: bash

   sheriff ${bam_} ${ref_} ${cells_} ${gtf_} \
       -cpu ${cpu} \
       --cnv_file ${cnv_} \
       --blacklist_file ${blacklist_} \
       --blacklist_seqs ${blacklist_seqs} \
       --edit_site_min_cells ${min_} \
       -o ${out_dir}

**Expected Output:**

.. code-block:: text

   Counting barcoded edits...
   Processed barcoded edits in 0.123 minutes

   Filtering canonical edits to those with criteria: min_cells: 1
   42 / 50 kept after min cells criteria

   Counting Non-barcoded edits...
   ...

Step 5: Examine Outputs
------------------------

Sheriff generates 13 output files in the output directory:

**Key Output Files:**

.. code-block:: bash

   ls ${out_dir}

   # Main results
   edit_site_info.txt                                    # Edit site metadata
   cell_allelic_dosage.canonical-edit-sites.parquet.gz   # Cell × edit-site counts
   cell_gene_mrna_counts.parquet.gz                      # Cell × gene expression

   # Supporting files
   edit_sites.bed                                        # BED format edit sites
   t7_barcode_edits.tsv                                  # Detailed edit events
   t7_barcoded_counts.parquet.gz                         # Barcoded UMI counts
   t7_nonbarcoded_counts.parquet.gz                      # Non-barcoded UMI counts

See :doc:`outputs` for detailed descriptions of all output files.

Step 6: Load Results in Python
-------------------------------

Read the parquet files using pandas:

.. code-block:: python

   import pandas as pd

   # Load allelic dosage (cell × edit-site)
   allelic_dosage = pd.read_parquet(
       "subset_500_cell_sheriff_output/cell_allelic_dosage.canonical-edit-sites.parquet.gz"
   )

   # Load gene expression (cell × gene)
   gene_counts = pd.read_parquet(
       "subset_500_cell_sheriff_output/cell_gene_mrna_counts.parquet.gz"
   )

   # Load edit site information
   edit_info = pd.read_csv(
       "subset_500_cell_sheriff_output/edit_site_info.txt",
       sep='\t'
   )

   # Explore the data
   print(f"Cells: {allelic_dosage.shape[0]}")
   print(f"Edit sites: {allelic_dosage.shape[1]}")
   print(f"Genes: {gene_counts.shape[1]}")

   # Get edit sites with > 10 cells
   edit_counts = (allelic_dosage > 0).sum(axis=0)
   major_edits = edit_counts[edit_counts > 10]
   print(f"Edit sites in > 10 cells: {len(major_edits)}")

Step 7: Load Results in R
--------------------------

Read the parquet files using the arrow package:

.. code-block:: r

   # Install arrow if needed
   install.packages("arrow")
   library(arrow)

   # Load allelic dosage
   allelic_dosage <- as.data.frame(
       read_parquet("subset_500_cell_sheriff_output/cell_allelic_dosage.canonical-edit-sites.parquet.gz")
   )

   # Set row names
   rownames(allelic_dosage) <- allelic_dosage[, '__index_level_0__']
   allelic_dosage <- allelic_dosage[, 1:(ncol(allelic_dosage)-1)]

   # Load gene expression
   gene_counts <- as.data.frame(
       read_parquet("subset_500_cell_sheriff_output/cell_gene_mrna_counts.parquet.gz")
   )
   rownames(gene_counts) <- gene_counts[, '__index_level_0__']
   gene_counts <- gene_counts[, 1:(ncol(gene_counts)-1)]

   # Explore
   cat(sprintf("Cells: %d\n", nrow(allelic_dosage)))
   cat(sprintf("Edit sites: %d\n", ncol(allelic_dosage)))
   cat(sprintf("Genes: %d\n", ncol(gene_counts)))

Next Steps
----------

* Explore :doc:`usage` for all parameter options
* Review :doc:`outputs` for detailed output descriptions
* Check :doc:`troubleshooting` if you encounter issues
* See the :doc:`api/modules` for programmatic usage

Important Notes
---------------

.. warning::
   The example BAM file is a filtered subset (200kb around edit sites). Edit site calling will be accurate, but gene expression counts will be incomplete due to filtered reads.

.. tip::
   For production analysis with full datasets, increase ``--cpu`` to speed up UMI counting and use ``--edit_site_min_cells 3`` as a quality filter.
