Changelog
=========

Version 1.1.3 (Current)
-----------------------

**Release Date:** February 2025

**Bug Fixes:**

* Fixed UMI counting bugs introduced during run-speed optimization in v1.1.0
* Corrected dtype in constructed arrays for mRNA UMI counting
* Improved error handling on threads

**Code Quality:**

* Added comprehensive requirements.txt
* Created pyproject.toml for modern packaging
* Fixed bare exception handlers
* Corrected typos in user-facing text
* Added package metadata to __init__.py
* Updated Python version requirement to >=3.10

**Documentation:**

* Added comprehensive Sphinx documentation
* Created CONTRIBUTING.md for contributors
* Improved code comments and docstrings

Version 1.1.1
-------------

**Release Date:** January 2025

**Changes:**

* Minor change in gene-counting logic to account for different split-pipe version outputs

Version 1.1.0
-------------

**Release Date:** December 2024

**Features:**

* Run-speed optimization
* Parallelization for gene counting
* Added multi-CPU support (``--cpu`` parameter)
* Implemented chunked genome processing

**Performance:**

* Significant speed improvements through Numba JIT compilation
* Reduced memory usage with sparse matrices
* Parallelized UMI counting across chromosomes

Version 1.0.0
-------------

**Release Date:** November 2024

**Initial Release:**

* Core functionality for CRISPR/Cas9 edit site detection
* Single-cell allelic dosage quantification
* Gene expression counting with T7 filtering
* K-mer based barcode matching
* Bidirectional edit site validation
* Support for blacklist/whitelist regions
* Copy number variant adjustment
* Multiple output formats (parquet, BAM, BED, TSV)

**Publication:**

Used for manuscript: "Joint single-cell profiling of CRISPR-Cas9 edits and transcriptomes reveals widespread off-target events and their effects on gene expression"

Roadmap
-------

Future enhancements under consideration:

**Version 1.2.x**

* Test suite with pytest
* CI/CD pipeline
* Logging module integration
* Pre-commit hooks

**Version 1.3.x**

* Type hints throughout codebase
* Improved error messages
* Extended documentation with more examples
* Tutorial notebooks

**Version 2.0.x**

* Support for additional CRISPR systems
* Multi-sample processing
* Integration with more single-cell pipelines
* Web-based quality control dashboard

Contributing
------------

See :doc:`contributing` for how to contribute to Sheriff development.

Migration Guides
----------------

1.0.0 → 1.1.x
^^^^^^^^^^^^^

**No breaking changes**. All command-line parameters remain compatible.

**Recommendations:**

* Use ``--cpu`` parameter for faster processing
* Verify UMI counts match expectations (bug fixes in 1.1.3)

1.1.x → 1.2.x (Future)
^^^^^^^^^^^^^^^^^^^^^^

Will be announced with release notes.
