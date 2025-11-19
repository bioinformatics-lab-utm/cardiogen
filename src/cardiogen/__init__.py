"""
Cardiogen - Cardiovascular Disease RNA-seq Analysis Pipeline
=============================================================

A space-efficient batch processing pipeline for RNA-seq data analysis.

Modules:
    sra_downloader: SRA data download and conversion utilities
    fastq_qc: Quality control and filtering modules
    
Example:
    >>> from cardiogen.sra_downloader import download_sra_files
    >>> from cardiogen.fastq_qc import QualityControl
"""

__version__ = "1.0.0"
__author__ = "Nicolae Drabcinski"
__email__ = "nicolae@example.com"

from pathlib import Path

# Package root directory
PACKAGE_ROOT = Path(__file__).parent.parent.parent
