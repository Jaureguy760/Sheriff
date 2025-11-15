About Sheriff
=============

Overview
--------

Sheriff is a bioinformatics tool developed at the Salk Institute for Biological Studies to analyze CRISPR/Cas9 editing outcomes in single cells using Superb-seq data.

The Project
-----------

**Goal:** Enable researchers to simultaneously measure CRISPR editing events and gene expression changes at single-cell resolution.

**Method:** Superb-seq combines CRISPR/Cas9 editing with single-cell RNA sequencing by detecting T7 barcode sequences inserted at edit sites.

**Output:** Comprehensive analysis of:

* Edit site locations and frequencies
* Allelic dosage (edited vs. unedited)
* Gene expression profiles
* Integration with standard scRNA-seq workflows

Scientific Background
---------------------

CRISPR/Cas9 technology enables precise genome editing, but understanding editing outcomes and their effects on gene expression requires high-throughput single-cell analysis.

**Superb-seq** (Single-cell profiling of CRISPR Edits with RNA-seq barcodes) addresses this by:

1. Incorporating T7 barcodes into donor DNA
2. Detecting inserted T7 sequences in scRNA-seq reads
3. Mapping edit sites genome-wide
4. Quantifying editing efficiency per cell
5. Measuring transcriptional consequences

Publication
-----------

Sheriff is described in:

   **Joint single-cell profiling of CRISPR-Cas9 edits and transcriptomes reveals widespread off-target events and their effects on gene expression**

   Michael H. Lorenzini, Brad Balderson, Karthyayani Sajeev, Aaron J. Ho, Graham McVicker

   *bioRxiv* 2025.02.07.636966

   doi: https://doi.org/10.1101/2025.02.07.636966

Key Findings
^^^^^^^^^^^^

* Comprehensive mapping of on-target and off-target CRISPR edits
* Single-cell resolution of editing outcomes
* Transcriptional impact of edit events
* Validation of Superb-seq methodology

The Team
--------

**Developers:**

* **Brad Balderson** - Lead Developer, Salk Institute

  Email: bbalderson@salk.edu

* **Michael Lorenzini** - Co-developer, Methodology

* **Aaron Ho** - Co-developer, Bioinformatics

**Principal Investigator:**

* **Graham McVicker** - Salk Institute for Biological Studies

Affiliations
------------

**Salk Institute for Biological Studies**

10010 North Torrey Pines Road

La Jolla, CA 92037, USA

https://www.salk.edu/

Funding
-------

This work was supported by:

* The Salk Institute for Biological Studies
* [Add specific grants if applicable]

Technology Stack
----------------

Sheriff is built with:

**Core:**

* Python 3.10+
* NumPy, Pandas, Polars for data manipulation
* Numba for JIT compilation

**Bioinformatics:**

* pysam for BAM file processing
* gtfparse for annotation parsing
* pyranges for genomic intervals
* BioPython for sequence analysis

**Performance:**

* FAISS for efficient k-mer matching
* SciPy sparse matrices for memory efficiency
* Multiprocessing for parallelization

**Packaging:**

* Setuptools, pip for distribution
* Sphinx for documentation
* ReadTheDocs theme for website

Community
---------

**GitHub:** https://github.com/BradBalderson/Sheriff

**Issues:** https://github.com/BradBalderson/Sheriff/issues

**Discussions:** [Coming soon]

**Contributors:** See CONTRIBUTING.md

License
-------

Sheriff is released under the BSD License.

See :doc:`license` for full text.

Citation
--------

If you use Sheriff in your research, please cite:

.. code-block:: bibtex

   @article{lorenzini2025sheriff,
     title={Joint single-cell profiling of CRISPR-Cas9 edits and transcriptomes reveals widespread off-target events and their effects on gene expression},
     author={Lorenzini, Michael H. and Balderson, Brad and Sajeev, Karthyayani and Ho, Aaron J. and McVicker, Graham},
     journal={bioRxiv},
     year={2025},
     doi={10.1101/2025.02.07.636966}
   }

Acknowledgments
---------------

We thank:

* The Salk Institute community for support
* Parse Biosciences for split-pipe software
* The open-source bioinformatics community
* All contributors and users

Contact
-------

For questions, bug reports, or collaboration inquiries:

* **Email:** bbalderson@salk.edu
* **GitHub Issues:** https://github.com/BradBalderson/Sheriff/issues
* **GitHub:** @BradBalderson

Related Projects
----------------

**Superb-seq Analysis:**

https://github.com/BradBalderson/superb_analysis

**Parse Biosciences (split-pipe):**

https://www.parsebiosciences.com/

**Scanpy (scRNA-seq analysis):**

https://scanpy.readthedocs.io/

**Seurat (R scRNA-seq):**

https://satijalab.org/seurat/
