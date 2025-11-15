"""
Sheriff - CRISPR/Cas9 Edit Site Detection in Single Cells

Sheriff processes aligned Superb-seq data to call edit sites and quantify
gene expression in single cells.
"""

__version__ = "1.1.3"
__author__ = "Brad Balderson, Michael Lorenzini, Aaron Ho"
__email__ = "bbalderson@salk.edu"
__license__ = "BSD"

# Import main functions for easier access
from .count_t7 import run_count_t7

__all__ = ['run_count_t7', '__version__']
