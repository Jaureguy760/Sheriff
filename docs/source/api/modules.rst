API Reference
=============

Sheriff can be used programmatically from Python scripts. This section documents the public API.

Main Modules
------------

.. toctree::
   :maxdepth: 2

   count_t7
   helpers
   rust

Quick Start
-----------

**Basic usage:**

.. code-block:: python

   from sheriff.count_t7 import run_count_t7

   # Run Sheriff programmatically
   run_count_t7(
       bam_file="aligned.bam",
       ref_file="reference.fa",
       barcode_file="barcodes.txt",
       gtf_file="genes.gtf",
       t7_barcode="GGGAGAGTAT",
       k=6,
       edit_dist=140,
       edit_site_min_cells=3,
       outdir="./results/",
       verbosity=1,
       n_cpus=4
   )

**Loading results:**

.. code-block:: python

   import pandas as pd

   # Load allelic dosage
   dosage = pd.read_parquet("results/cell_allelic_dosage.canonical-edit-sites.parquet.gz")

   # Load gene expression
   expression = pd.read_parquet("results/cell_gene_mrna_counts.parquet.gz")

   # Load edit info
   edit_info = pd.read_csv("results/edit_site_info.txt", sep='\t')

Package Information
-------------------

.. code-block:: python

   import sheriff

   print(sheriff.__version__)  # '1.1.3'
   print(sheriff.__author__)   # 'Brad Balderson, Michael Lorenzini, Aaron Ho'

See Also
--------

* :doc:`../tutorial` - Step-by-step examples
* :doc:`../outputs` - Output file descriptions
