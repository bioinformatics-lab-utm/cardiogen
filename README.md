# 🧬 Cardiogen

[![Python](https://img.shields.io/badge/Python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> Space-efficient batch processing pipeline for cardiovascular disease RNA-seq data analysis

## 🎯 Features

- 🚀 **Batch Processing** - Process 1000+ samples with minimal disk space
- 💾 **Space Efficient** - 99% disk space savings through automatic cleanup
- 📊 **Quality Control** - Comprehensive FastQC analysis for all samples
- 📈 **Detailed Reports** - JSON, CSV, and HTML reports with all metrics
- 🔄 **Resumable** - Continue processing from any interrupted batch
- 🏥 **Disease-specific** - Organized results by cardiovascular disease type
- 🔬 **Metadata Download** - Fetch all patient/clinical metadata from NCBI SRA

## 📋 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/cardiogen.git
cd cardiogen

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install system tools
sudo apt-get install sra-toolkit fastqc  # Ubuntu/Debian
```

### Basic Usage

```bash
# Test run (2 samples)
python scripts/run_pipeline.py --max-samples 2 --batch-size 1

# Process specific disease
python scripts/run_pipeline.py --disease "Aortic Valve Disease" --batch-size 5

# Process all data
python scripts/run_pipeline.py --batch-size 10
```

### Download Metadata (Patient/Clinical Data)

Download all available metadata including patient sex, age, disease status, tissue type, etc:

```bash
# Download all metadata for all SRR IDs in CSV files
python scripts/download_all_metadata.py --email your_email@example.com

# Run in background
nohup python scripts/download_all_metadata.py --email your_email@example.com > metadata_download.log 2>&1 &

# Monitor progress
./scripts/monitor_metadata.sh
```

**Output**: Complete CSV file with 60+ columns including all patient/clinical attributes.  
**See**: [Metadata Download Guide](docs/METADATA_DOWNLOAD.md) for details.

## 📁 Project Structure

```
cardiogen/
├── src/cardiogen/           # Core modules
│   ├── __init__.py
│   ├── sra_downloader.py    # SRA download utilities
│   └── fastq_qc.py          # Quality control module
│
├── scripts/                 # Executable scripts
│   ├── run_pipeline.py      # Main pipeline script
│   ├── eda_analysis.py      # Exploratory data analysis
│   ├── download_all_metadata.py  # Download SRA metadata
│   ├── monitor_pipeline.sh  # Monitor pipeline progress
│   ├── monitor_metadata.sh  # Monitor metadata download
│   └── cleanup.sh           # Cleanup utility
│
├── config/                  # Configuration files
│   └── config.py            # Pipeline configuration
│
├── data/                    # Data directory (gitignored)
│   ├── qc/                 # FastQC reports
│   ├── results/            # Summary results
│   └── *.csv               # Input datasets
│
├── docs/                    # Documentation
│   ├── GUIDE.md
│   └── README_BATCH.md
│
├── tests/                   # Unit tests
├── examples/                # Usage examples
│
├── README.md
├── requirements.txt
├── setup.py
└── LICENSE
```

## 💾 Disk Space Efficiency

**Traditional Approach:**
- 100 samples × 5GB = 500GB required

**Cardiogen Batch Pipeline:**
- Max simultaneous: 10 samples × 5GB = 50GB
- Final storage: 100 samples × 50MB QC = 5GB
- **Space savings: 99%!** 🎉

## 📊 Supported Diseases

- Aortic Valve Disease
- Marfan Syndrome
- Loeys-Dietz Syndrome
- Short QT Syndrome
- Familial Hypercholesterolemia

## 📖 Documentation

- [Installation Guide](docs/GUIDE.md#installation)
- [Usage Examples](EXAMPLES.md)
- [Metadata Download Guide](docs/METADATA_DOWNLOAD.md) - Download all patient/clinical data
- [API Documentation](docs/GUIDE.md#api)
- [Contributing Guidelines](CONTRIBUTING.md)

## 🔧 Configuration

Edit `config/config.py` to customize pipeline behavior:

```python
# Batch processing settings
MAX_PARALLEL_PROCESSES = 4
BATCH_SIZE = 10

# Quality thresholds
QUALITY_THRESHOLD = 20
LENGTH_THRESHOLD = 50
GC_MIN = 30
GC_MAX = 70
```

## 📈 Example Analysis

```python
import pandas as pd

# Load results
df = pd.read_csv('data/results/processing_summary.csv')

# Analyze by disease
disease_stats = df.groupby('disease')['total_sequences'].describe()
print(disease_stats)

# Plot GC content distribution
import matplotlib.pyplot as plt
df.boxplot(column='gc_content', by='disease', figsize=(12,6))
plt.savefig('gc_content_analysis.png')
```

## 🧪 Testing

```bash
# Run tests
pytest tests/

# With coverage
pytest --cov=src tests/
```

## 🤝 Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 👥 Authors

- Nicolae Drabcinski

## 🙏 Acknowledgments

- [SRA Toolkit](https://github.com/ncbi/sra-tools) by NCBI
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) by Babraham Bioinformatics
- [MultiQC](https://multiqc.info/) by Phil Ewels

## 📧 Contact

For questions or issues, please [open an issue](https://github.com/YOUR_USERNAME/cardiogen/issues) on GitHub.

---

**Made with ❤️ for cardiovascular disease research**

## 📋 Overview

This pipeline processes large-scale SRA (Sequence Read Archive) datasets efficiently by:
- **Batch processing** - Processes samples in small batches to minimize disk usage
- **Automatic cleanup** - Removes raw FASTQ files after QC analysis
- **Quality control** - Runs FastQC on all samples
- **Space-efficient** - Saves 99% disk space by keeping only QC reports

## 🎯 Key Features

- ✅ **Space-efficient batch processing** - Process 1000+ samples with only 50-100GB free space
- ✅ **Automatic data management** - Downloads, processes, and cleans up automatically
- ✅ **Comprehensive QC** - FastQC analysis for all samples
- ✅ **Detailed reports** - JSON, CSV, and text reports with all metrics
- ✅ **Disease-based organization** - Results organized by cardiovascular disease type
- ✅ **Resumable** - Can resume processing from any batch if interrupted

## 📊 Supported Diseases

- Aortic Valve Disease
- Marfan Syndrome
- Loeys-Dietz Syndrome
- Short QT Syndrome
- Familial Hypercholesterolemia

## 🚀 Quick Start

### Prerequisites

**Required tools:**
- Python 3.11+
- SRA Toolkit (prefetch, fasterq-dump)
- FastQC
- MultiQC (optional)

**Installation:**
```bash
# Install SRA Toolkit
sudo apt-get install sra-toolkit

# Install FastQC
sudo apt-get install fastqc

# Install Python dependencies
python -m venv .venv
source .venv/bin/activate
pip install pandas biopython multiqc
```

### Basic Usage

```bash
# Activate virtual environment
source .venv/bin/activate

# Test run (2 samples)
python batch_pipeline.py --max-samples 2 --batch-size 1

# Process specific disease
python batch_pipeline.py --disease "Aortic Valve Disease" --batch-size 5

# Process all data
python batch_pipeline.py --batch-size 10
```

### Background Processing

```bash
# Run in background (recommended for large datasets)
nohup python batch_pipeline.py --batch-size 10 > output.log 2>&1 &

# Monitor progress
tail -f output.log
```

## 📁 Project Structure

```
cardiogen/
├── batch_pipeline.py          # Main batch processing script ⭐
├── pipeline_config.py         # Configuration settings
├── sra_downloader.py          # SRA download module
├── fastq_qc.py               # Quality control module
├── sorted_sra_samples_by_disease_MOROZ_file.csv  # Input data
│
├── data/                      # Processing output
│   ├── qc/                   # FastQC reports (KEPT)
│   │   ├── aortic_valve_disease/
│   │   ├── marfan_syndrome/
│   │   └── ...
│   ├── results/              # Summary results (KEPT)
│   │   ├── batch_*_results.json
│   │   ├── processing_summary.csv
│   │   └── summary_report.txt
│   └── temp_batch_*/         # Temporary files (AUTO-DELETED)
│
├── docs/                      # Documentation
│   ├── GUIDE.md              # Complete guide
│   └── README_BATCH.md       # Batch processing details
│
├── old_scripts/              # Legacy scripts
└── archive/                  # Archived data
```

## 💾 Disk Space Management

**Traditional approach:**
- 100 samples × 5GB = **500GB required** 

**Batch pipeline:**
- Processes 10 samples at a time = **50GB max**
- Keeps only QC reports = **5GB final**
- **99% space savings!** 🎉

## 📈 Processing Time

**Per sample:** ~2-3 minutes
- Download: 30-60 seconds
- Conversion: 30-90 seconds  
- FastQC: 30-60 seconds

**Full dataset (2394 samples):**
- With `--batch-size 10`: ~8-12 hours
- With `--batch-size 5`: ~10-15 hours

## 📊 Output Files

### QC Reports
- `data/qc/{disease}/{accession}/` - FastQC HTML reports and data
- Visual quality assessment for each sample

### Summary Files
- `processing_summary.csv` - Tabular data with all metrics
- `summary_report.txt` - Human-readable summary
- `batch_*_results.json` - Detailed batch results

### Example Analysis

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('data/results/processing_summary.csv')

# Analyze by disease
print(df.groupby('disease')['total_sequences'].describe())

# Plot GC content distribution
df.boxplot(column='gc_content', by='disease', figsize=(12,6))
plt.savefig('gc_by_disease.png')
```

## 🔧 Configuration

Edit `pipeline_config.py` to customize:

```python
# Batch size (samples processed simultaneously)
MAX_PARALLEL_PROCESSES = 4

# Quality thresholds
QUALITY_THRESHOLD = 20
LENGTH_THRESHOLD = 50

# Keep raw files (not recommended)
keep_raw_fastq = False  # Set to True to keep FASTQ files
```

## 🐛 Troubleshooting

### Common Issues

**"No space left on device"**
```bash
# Clean temporary files
rm -rf data/temp_batch_*

# Reduce batch size
python batch_pipeline.py --batch-size 2
```

**"Command not found: python"**
```bash
source .venv/bin/activate
```

**Missing tools**
```bash
# Check installed tools
which prefetch fasterq-dump fastqc

# Install missing tools
sudo apt-get install sra-toolkit fastqc
```

## 📚 Documentation

- [Complete Guide](docs/GUIDE.md) - Detailed documentation
- [Batch Processing](docs/README_BATCH.md) - Batch pipeline specifics
- [Configuration](pipeline_config.py) - Configuration options

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## 👥 Authors

- Nicolae Drabcinski

## 🙏 Acknowledgments

- SRA Toolkit by NCBI
- FastQC by Babraham Bioinformatics
- MultiQC by Phil Ewels

## 📧 Contact

For questions or issues, please open an issue on GitHub.

---

**Note:** This pipeline is designed for research purposes. Always verify results with domain experts.
