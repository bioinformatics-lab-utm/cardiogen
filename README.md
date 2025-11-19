# 🧬 CardioGen - Cardiovascular Disease RNA-seq Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-bioinformatics--lab--utm-lightgrey)](https://github.com/bioinformatics-lab-utm/cardiogen)

> Comprehensive SRA metadata analysis and batch processing pipeline for cardiovascular disease genomics research

## 📊 Dataset Overview

- **Total Samples:** 4,596 SRA runs
- **Illumina Filtered:** 2,880 runs from 2,712 unique patients
- **Metadata Columns:** 300+ attributes per sample
- **Phenotype Data:** 27 clinical/demographic columns analyzed
- **Visualizations:** 45+ publication-quality graphs
- **Storage:** 5.1TB processed data (3.3TB SRA, 1.8TB archive, 650MB QC)

## 🎯 Key Features

### Data Analysis
- 📈 **Comprehensive Metadata Analysis** - Complete exploration of 300+ metadata columns
- 👥 **Sample-level Deduplication** - Patient-level analysis (2,712 unique from 2,880 runs)
- 🧬 **Phenotype Analysis** - 10 detailed graphs covering sex, age, disease, tissue, treatment, conditions
- 📊 **45+ Visualizations** - Publication-ready figures (300 DPI, multiple formats)
- 📑 **PowerPoint Generator** - Automated presentation creation (50 slides)

### Data Processing
- 🚀 **Batch Pipeline** - Process 1000+ samples efficiently
- 💾 **Space Optimization** - Automatic cleanup after processing
- 📊 **Quality Control** - FastQC analysis for all samples
- 🔄 **Resumable Processing** - Continue from interrupted batches
- 🏥 **Disease Organization** - Results organized by cardiovascular disease type

### Technical
- 🐍 **Python 3.11+** - Modern OOP design with type hints
- 📓 **Jupyter Integration** - Interactive analysis notebooks
- 🔧 **Modular Architecture** - Reusable components in `src/cardiogen/`
- 📦 **Full Dependency Management** - Virtual environment with requirements.txt

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/bioinformatics-lab-utm/cardiogen.git
cd cardiogen

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install system tools (Ubuntu/Debian)
sudo apt-get install sra-toolkit fastqc
```

### Basic Usage

#### Metadata Analysis
```bash
# Run complete metadata analysis
python sra_metadata_analysis.py

# Generate PowerPoint presentation
python create_presentation.py
```

#### Data Download Pipeline
```bash
# Test run (2 samples)
python scripts/run_pipeline.py --max-samples 2 --batch-size 1

# Process specific disease
python scripts/run_pipeline.py --disease "Aortic Valve Disease" --batch-size 5

# Process all data with batch processing
python scripts/run_pipeline.py --batch-size 10
```

#### View Metadata
```bash
# Download and view SRA metadata
python scripts/download_all_metadata.py

# View specific samples
python scripts/view_metadata.py

# Interactive EDA
python scripts/metadata_eda.py
```

## 📁 Project Structure

```
cardiogen/
├── sra_metadata_analysis.py    # Main analysis script (103KB, 45+ graphs)
├── create_presentation.py       # PowerPoint generator (50 slides)
├── requirements.txt             # Python dependencies
├── config/
│   └── config.py               # Configuration settings
├── src/cardiogen/              # Core modules
│   ├── sra_downloader.py       # SRA data download
│   └── fastq_qc.py             # Quality control
├── scripts/                     # Utility scripts
│   ├── run_pipeline.py         # Batch processing pipeline
│   ├── download_all_metadata.py
│   ├── metadata_eda.py
│   └── *.sh                    # Monitoring scripts
├── data/
│   ├── metadata/               # SRA metadata (CSV/JSON)
│   ├── analysis_results/       # Results and figures
│   │   └── figures/            # 45+ PNG visualizations
│   ├── adapters/               # Trimmomatic adapters
│   ├── qc/                     # FastQC reports
│   └── results/                # Processing results
├── SRA_Complete_Phenotypes.pptx # Generated presentation (8.2MB)
└── PROJECT_AUDIT.md            # Comprehensive project audit
```

## 📊 Analysis Outputs

### Visualizations (45 graphs)
1. **Instruments & Technology**
   - Sequencing platform distribution
   - Library strategy/layout/selection
   - Sequencing depth analysis

2. **Demographics & Phenotypes**
   - Sex distribution (detailed & aggregated)
   - Age distribution histogram
   - Metadata completeness heatmaps

3. **Clinical Data**
   - Disease distribution (10+ categories)
   - Tissue/body site analysis
   - Treatment and condition data
   - Disease state tracking

4. **Sample-level Analysis**
   - Runs per sample distribution
   - Patient deduplication (2,712 unique)
   - Top samples by technical replicates

5. **Quality Metrics**
   - Data completeness by category
   - Statistical summary heatmaps
   - Correlation matrices

### Generated Files
- **CSV Reports:** Summary statistics, filtered metadata, quality metrics
- **PNG Figures:** 300 DPI publication-ready visualizations
- **PowerPoint:** Automated presentation with all graphs
- **JSON Metadata:** Complete SRA records with full hierarchy

## 🧬 Metadata Analysis

### SRA Hierarchy
```
STUDY → SAMPLE → EXPERIMENT → RUN
```
- **Study:** Research project (e.g., PRJNA123456)
- **Sample:** Biological specimen (patient)
- **Experiment:** Library preparation
- **Run:** Sequencing run (technical replicate)

### Key Insights
- **2,712 unique patients** from 2,880 Illumina runs
- **87 patients** with multiple technical replicates
- **27 phenotype columns** analyzed (sex, age, disease, tissue, treatment, etc.)
- **300+ total metadata attributes** per sample

## 🔧 Advanced Usage

### Custom Analysis
```python
from sra_metadata_analysis import SRAMetadataAnalyzer

# Initialize analyzer
analyzer = SRAMetadataAnalyzer('data/metadata/sra_metadata_complete.csv')

# Run specific analyses
analyzer.analyze_comprehensive_phenotypes()
analyzer.analyze_sample_level_summary()
analyzer.analyze_instruments()
```

### Presentation Customization
```python
from create_presentation import SRAPresentationCreator

# Create custom presentation
creator = SRAPresentationCreator()
creator.create_presentation()  # Generates SRA_Complete_Phenotypes.pptx
```

### Pipeline Configuration
Edit `config/config.py` to customize:
- Data directories
- Batch sizes
- Quality control thresholds
- Output formats

## 📈 Performance

- **Analysis Runtime:** ~5-10 minutes for full metadata analysis
- **Visualization Generation:** ~2-3 minutes for 45 graphs
- **Presentation Creation:** ~30 seconds
- **Data Processing:** ~10-15 minutes per batch (depending on batch size)

## 🐛 Troubleshooting

### Common Issues

**Missing Dependencies**
```bash
pip install pandas numpy matplotlib seaborn scipy python-pptx
```

**SRA Tools Not Found**
```bash
# Ubuntu/Debian
sudo apt-get install sra-toolkit

# MacOS
brew install sra-tools
```

**Disk Space Issues**
- Check available space: `df -h`
- Review archive directory: `du -sh archive/`
- Clean temporary files: `bash scripts/cleanup.sh`

## 📝 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{cardiogen2025,
  title = {CardioGen: Cardiovascular Disease RNA-seq Analysis Pipeline},
  author = {Bioinformatics Lab UTM},
  year = {2025},
  url = {https://github.com/bioinformatics-lab-utm/cardiogen}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📧 Contact

**Bioinformatics Lab UTM**
- GitHub: [@bioinformatics-lab-utm](https://github.com/bioinformatics-lab-utm)
- Repository: [cardiogen](https://github.com/bioinformatics-lab-utm/cardiogen)

## 🙏 Acknowledgments

- NCBI SRA for cardiovascular disease datasets
- FastQC for quality control tools
- Python scientific computing community (pandas, matplotlib, seaborn)

---

**Last Updated:** November 19, 2025