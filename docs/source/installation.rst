Installation
============

Prerequisites
-------------

* Python 3.10 or higher
* Conda or Mamba (recommended for dependency management)
* Git

Recommended Installation
------------------------

Using Mamba (Faster)
^^^^^^^^^^^^^^^^^^^^

Mamba is a faster alternative to Conda. If you don't have it, install from `mamba-org <https://mamba.readthedocs.io/>`_.

.. code-block:: bash

   # Create a new environment
   mamba create -n sheriff_env python=3.10
   mamba activate sheriff_env

   # Install dependencies
   mamba install conda-forge::scipy conda-forge::numpy==1.26.4 \
       bioconda::gtfparse conda-forge::faiss-cpu conda-forge::numba \
       conda-forge::biopython=1.81 typing_extensions typer \
       bioconda::pyranges bioconda::pysam conda-forge::polars

   # Clone and install Sheriff
   git clone https://github.com/BradBalderson/Sheriff.git
   cd Sheriff
   pip install .

   # Verify installation
   sheriff --help

Using Conda
^^^^^^^^^^^

If you prefer Conda:

.. code-block:: bash

   # Replace 'mamba' with 'conda' in the commands above
   conda create -n sheriff_env python=3.10
   conda activate sheriff_env
   # ... (same as above but with 'conda' instead of 'mamba')

From PyPI (Coming Soon)
^^^^^^^^^^^^^^^^^^^^^^^

.. note::
   PyPI installation is not yet available. Use the Git installation method above.

Development Installation
------------------------

For contributing to Sheriff:

.. code-block:: bash

   # Clone your fork
   git clone https://github.com/YOUR_USERNAME/Sheriff.git
   cd Sheriff

   # Create environment
   mamba create -n sheriff_dev python=3.10
   mamba activate sheriff_dev

   # Install dependencies
   pip install -r requirements.txt

   # Install development dependencies
   pip install pytest pytest-cov black ruff mypy sphinx sphinx-rtd-theme

   # Install in editable mode
   pip install -e .

Verify Installation
-------------------

Check that Sheriff is correctly installed:

.. code-block:: bash

   sheriff --version
   sheriff --help

System-Specific Notes
---------------------

macOS
^^^^^

On macOS, you may need to set up library paths:

.. code-block:: bash

   export LDFLAGS="-L$CONDA_PREFIX/lib"
   export CPPFLAGS="-I$CONDA_PREFIX/include"
   export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

Add these to your ``~/.bashrc`` or ``~/.zshrc`` to make them permanent.

Linux
^^^^^

Sheriff has been tested on Rocky 9 and Ubuntu distributions. No special configuration needed.

Windows
^^^^^^^

Sheriff has not been tested on Windows. We recommend using WSL (Windows Subsystem for Linux).

Troubleshooting
---------------

zlib Not Found
^^^^^^^^^^^^^^

If you encounter zlib errors:

.. code-block:: bash

   mamba install conda-forge::zlib
   export LDFLAGS="-L$CONDA_PREFIX/lib"
   export CPPFLAGS="-I$CONDA_PREFIX/include"

Numba Compatibility
^^^^^^^^^^^^^^^^^^^

Sheriff requires ``numba>=0.56.0``. If you have issues:

.. code-block:: bash

   mamba install conda-forge::numba --force-reinstall

FAISS Installation Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^

For FAISS CPU issues:

.. code-block:: bash

   mamba install conda-forge::faiss-cpu=1.10.0 --force-reinstall

Getting Help
------------

If you encounter installation issues:

1. Check the :doc:`troubleshooting` guide
2. Search `existing issues <https://github.com/BradBalderson/Sheriff/issues>`_
3. Open a new issue with your error message and system information
4. Email: bbalderson@salk.edu
