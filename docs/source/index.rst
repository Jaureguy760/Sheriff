Sheriff Documentation
=====================

.. image:: https://img.shields.io/badge/python-3.10+-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python Version

.. image:: https://img.shields.io/badge/License-BSD-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: License

**Sheriff** is a bioinformatics tool for identifying CRISPR/Cas9 edit sites in single cells and quantifying gene expression from Superb-seq data.

Sheriff processes aligned Superb-seq data to:

* Call CRISPR/Cas9 edit sites with single-cell resolution
* Quantify allelic dosage (edited vs. unedited alleles)
* Generate UMI-based gene expression counts
* Filter T7-related reads to improve expression accuracy

Key Features
------------

**End-to-End Analysis**
   Process BAM files from split-pipe to generate edit site calls and gene expression matrices

**High Performance**
   Optimized with Numba JIT compilation, optional Rust acceleration (10-50x faster BAM filtering), and parallel processing for fast analysis

**Rust Acceleration (NEW!)**
   Optional Rust implementation for BAM filtering provides 10-50x speedup on performance-critical operations with automatic fallback to Python

**Comprehensive Output**
   13 different output files including edit sites, allelic dosage, and gene expression matrices

**Flexible Parameters**
   Customizable filtering criteria for edit site calling and quality control

**Professional Documentation**
   Complete API reference with auto-generated documentation from code docstrings

Quick Start
-----------

Installation
^^^^^^^^^^^^

.. code-block:: bash

   # Create conda environment
   mamba create -n sheriff_env python=3.10
   mamba activate sheriff_env

   # Install dependencies
   mamba install conda-forge::scipy conda-forge::numpy==1.26.4 \
       bioconda::gtfparse conda-forge::faiss-cpu conda-forge::numba \
       conda-forge::biopython=1.81 typing_extensions typer \
       bioconda::pyranges bioconda::pysam

   # Install Sheriff
   git clone https://github.com/BradBalderson/Sheriff.git
   cd Sheriff
   pip install .

Basic Usage
^^^^^^^^^^^

.. code-block:: bash

   sheriff <bam_file> <ref_file> <barcode_file> <gtf_file> \
       --cpu 4 \
       --cnv_file <cnv_file> \
       --blacklist_file <blacklist_file> \
       --edit_site_min_cells 3 \
       -o <output_directory>

See the :doc:`tutorial` for a complete example with real data.

Citation
--------

If you use Sheriff in your research, please cite:

   Joint single-cell profiling of CRISPR-Cas9 edits and transcriptomes reveals widespread off-target events and their effects on gene expression.
   Michael H. Lorenzini, Brad Balderson, Karthyayani Sajeev, Aaron J. Ho, Graham McVicker.
   *bioRxiv* 2025.02.07.636966; doi: https://doi.org/10.1101/2025.02.07.636966

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   tutorial
   usage
   outputs
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/modules
   api/count_t7
   api/helpers

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   changelog

.. toctree::
   :maxdepth: 1
   :caption: About

   about
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
