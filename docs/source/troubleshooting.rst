Troubleshooting
===============

Common issues and solutions for Sheriff.

Installation Issues
-------------------

zlib Not Found
^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   OSError: libz.so.1: cannot open shared object file

**Solution:**

.. code-block:: bash

   mamba install conda-forge::zlib
   export LDFLAGS="-L$CONDA_PREFIX/lib"
   export CPPFLAGS="-I$CONDA_PREFIX/include"
   export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

Add the exports to your ``~/.bashrc`` for persistence.

Numba Compilation Errors
^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   NumbaTypeSafetyWarning: ...

**Solution:**

.. code-block:: bash

   mamba install conda-forge::numba=0.56.4 --force-reinstall

FAISS Not Found
^^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   ModuleNotFoundError: No module named 'faiss'

**Solution:**

.. code-block:: bash

   mamba install conda-forge::faiss-cpu=1.10.0

Runtime Errors
--------------

BAM File Not Found or Invalid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   ValueError: could not open BAM file

**Checklist:**

1. Verify BAM file exists and is readable
2. Check BAM is sorted: ``samtools view -H file.bam | grep SO:``
3. Index the BAM: ``samtools index file.bam``
4. Verify BAM has required tags (CB, UB):

.. code-block:: bash

   samtools view file.bam | head -1 | grep -o "CB:Z:[^ ]*"

Reference Genome Mismatch
^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   Warning could not load from reference sequence

**Solutions:**

1. Ensure reference matches alignment genome
2. Check chromosome naming (chr1 vs 1)
3. Verify FASTA is indexed:

.. code-block:: bash

   samtools faidx reference.fa
   ls reference.fa.fai  # Should exist

Cell Barcode Issues
^^^^^^^^^^^^^^^^^^^

**Problem:** No cells detected or very few edit sites called.

**Solutions:**

1. Verify barcode format matches BAM:

.. code-block:: python

   import pysam
   bam = pysam.AlignmentFile("file.bam")
   read = next(bam)
   print(read.get_tag("CB"))  # Check format

2. Ensure barcode file has one barcode per line (no header)
3. Check for whitespace: ``cat barcodes.txt | wc -l``

Memory Issues
^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   MemoryError: Unable to allocate array

**Solutions:**

1. Reduce chunk size:

.. code-block:: bash

   sheriff ... --chunk 100  # Default is 250

2. Reduce number of CPUs (reduces parallel memory usage):

.. code-block:: bash

   sheriff ... --cpu 1

3. Process smaller cell subsets
4. Use a machine with more RAM

Performance Issues
------------------

Very Slow Processing
^^^^^^^^^^^^^^^^^^^^

**Solutions:**

1. Increase CPUs for parallelization:

.. code-block:: bash

   sheriff ... --cpu 8

2. Use blacklist to filter problematic regions:

.. code-block:: bash

   sheriff ... --blacklist blacklist.bed

3. Check BAM file is properly indexed
4. Verify disk I/O isn't bottleneck (use local SSD if possible)

High Memory Usage
^^^^^^^^^^^^^^^^^

**Solutions:**

1. Lower chunk size: ``--chunk 100``
2. Reduce CPU count: ``--cpu 2``
3. Process cells in batches

Output Issues
-------------

Empty Output Files
^^^^^^^^^^^^^^^^^^

**Problem:** No edit sites called.

**Checklist:**

1. Check verbosity output for warnings:

.. code-block:: bash

   sheriff ... -v 2 2>&1 | tee sheriff.log

2. Lower minimum cell threshold:

.. code-block:: bash

   sheriff ... --edit_site_min_cells 1

3. Disable bidirectional requirement:

.. code-block:: bash

   sheriff ... --no-bidirectional_inserts

4. Verify T7 barcode is correct: ``--t7 GGGAGAGTAT``

Parquet Files Won't Load
^^^^^^^^^^^^^^^^^^^^^^^^^

**Error:**

.. code-block:: text

   ArrowInvalid: Not a Parquet file

**Solutions:**

.. code-block:: bash

   # Verify file integrity
   file output.parquet.gz

   # Try rerunning Sheriff
   rm -rf output_directory/
   sheriff ... -o output_directory/

**Python loading:**

.. code-block:: python

   import pandas as pd
   df = pd.read_parquet("file.parquet.gz")  # .gz extension should work

**R loading:**

.. code-block:: r

   library(arrow)
   df <- read_parquet("file.parquet.gz")

Analysis Issues
---------------

Too Many Edit Sites Called
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Solutions:**

1. Increase minimum cells:

.. code-block:: bash

   sheriff ... --edit_site_min_cells 5

2. Use blacklist to filter known problematic regions
3. Enable/verify bidirectional requirement (default)
4. Reduce stranded edit distance:

.. code-block:: bash

   sheriff ... --stranded_edit_dist 10

Too Few Edit Sites Called
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Solutions:**

1. Lower minimum cells (exploration mode):

.. code-block:: bash

   sheriff ... --edit_site_min_cells 1

2. Disable bidirectional requirement:

.. code-block:: bash

   sheriff ... --no-bidirectional_inserts

3. Check if known sites are in whitelist
4. Increase verbosity to see filtering stats:

.. code-block:: bash

   sheriff ... -v 2

Debugging Tips
--------------

Enable Debug Logging
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   sheriff ... -v 2 2>&1 | tee debug.log

Check Intermediate Files
^^^^^^^^^^^^^^^^^^^^^^^^^

Examine the BAM subset files to verify read filtering:

.. code-block:: bash

   samtools view output/t7_barcoded_only.bam | less
   samtools flagstat output/t7_only.bam

Inspect Edit Details
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   import pandas as pd

   # Load detailed edit information
   edits = pd.read_csv("output/t7_barcode_edits.tsv", sep='\t')

   # Check edit distribution
   print(edits.groupby('edit_site').size().sort_values(ascending=False))

   # Inspect specific edit
   site = "chr1:12345"
   site_edits = edits[edits['edit_site'] == site]
   print(site_edits)

Getting Help
------------

If you can't resolve your issue:

1. **Check documentation:**

   * :doc:`installation`
   * :doc:`usage`
   * :doc:`outputs`

2. **Search existing issues:**

   https://github.com/BradBalderson/Sheriff/issues

3. **Open a new issue:**

   Include:

   * Sheriff version
   * Full command used
   * Complete error message
   * System information (OS, Python version)
   * Relevant log output (``-v 2``)

4. **Contact:**

   bbalderson@salk.edu
