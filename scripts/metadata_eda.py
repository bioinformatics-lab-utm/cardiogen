#!/usr/bin/env python3
"""
Exploratory Data Analysis for SRA Metadata.
Analyzes patient/clinical data from downloaded metadata.
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from collections import Counter

# Configuration
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10
sns.set_style('whitegrid')


class MetadataEDA:
    """EDA for SRA metadata."""
    
    def __init__(self, csv_file: str, output_dir: str = None):
        """
        Initialize EDA.
        
        Args:
            csv_file: Path to metadata CSV file
            output_dir: Directory to save reports and plots
        """
        self.csv_file = Path(csv_file)
        self.output_dir = Path(output_dir) if output_dir else self.csv_file.parent / 'eda'
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"Loading metadata from {self.csv_file}...")
        self.df = pd.read_csv(csv_file)
        print(f"Loaded {len(self.df)} records with {len(self.df.columns)} columns")
        
        # Categorize columns
        self.sample_cols = [c for c in self.df.columns if c.startswith('sample_')]
        self.run_cols = [c for c in self.df.columns if c.startswith('run_')]
        self.exp_cols = [c for c in self.df.columns if c.startswith(('experiment_', 'exp_attr_'))]
        self.study_cols = [c for c in self.df.columns if c.startswith('study_')]
        
    def generate_report(self):
        """Generate comprehensive EDA report."""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = self.output_dir / f'metadata_eda_report_{timestamp}.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("SRA METADATA - EXPLORATORY DATA ANALYSIS REPORT\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Source file: {self.csv_file.name}\n")
            f.write(f"Total samples: {len(self.df)}\n")
            f.write(f"Total columns: {len(self.df.columns)}\n")
            f.write("="*80 + "\n\n")
            
            # Dataset overview
            f.write("DATASET OVERVIEW\n")
            f.write("-"*80 + "\n")
            f.write(f"Sample/Patient attributes: {len(self.sample_cols)}\n")
            f.write(f"Run attributes: {len(self.run_cols)}\n")
            f.write(f"Experiment attributes: {len(self.exp_cols)}\n")
            f.write(f"Study attributes: {len(self.study_cols)}\n\n")
            
            # Demographics
            f.write("DEMOGRAPHICS\n")
            f.write("-"*80 + "\n")
            
            # Sex/Gender
            sex_cols = [c for c in self.sample_cols if any(x in c.lower() for x in ['sex', 'gender'])]
            for col in sex_cols:
                if col in self.df.columns:
                    f.write(f"\n{col}:\n")
                    counts = self.df[col].value_counts()
                    for val, count in counts.items():
                        f.write(f"  {val}: {count} ({count/len(self.df)*100:.1f}%)\n")
                    missing = self.df[col].isna().sum()
                    if missing > 0:
                        f.write(f"  Missing: {missing} ({missing/len(self.df)*100:.1f}%)\n")
            
            # Age
            age_cols = [c for c in self.sample_cols if 'age' in c.lower()]
            if age_cols:
                f.write(f"\nAge-related attributes ({len(age_cols)}):\n")
                for col in age_cols:
                    non_null = self.df[col].dropna()
                    if len(non_null) > 0:
                        f.write(f"  {col}: {len(non_null)} samples\n")
                        # Try to get statistics if numeric
                        try:
                            numeric_vals = pd.to_numeric(non_null, errors='coerce').dropna()
                            if len(numeric_vals) > 0:
                                f.write(f"    Mean: {numeric_vals.mean():.1f}, Median: {numeric_vals.median():.1f}\n")
                                f.write(f"    Range: {numeric_vals.min():.1f} - {numeric_vals.max():.1f}\n")
                        except:
                            # Show examples if not numeric
                            f.write(f"    Examples: {', '.join(map(str, non_null.head(3).values))}\n")
            
            # Disease/Phenotype
            f.write("\n\nDISEASE/PHENOTYPE\n")
            f.write("-"*80 + "\n")
            disease_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                           ['disease', 'phenotype', 'condition', 'diagnosis', 'affected'])]
            
            for col in disease_cols[:5]:
                if col in self.df.columns:
                    non_null = self.df[col].dropna()
                    if len(non_null) > 0:
                        f.write(f"\n{col}:\n")
                        counts = self.df[col].value_counts().head(10)
                        for val, count in counts.items():
                            val_str = str(val)[:60]
                            f.write(f"  {val_str}: {count}\n")
                        if len(self.df[col].value_counts()) > 10:
                            f.write(f"  ... and {len(self.df[col].value_counts()) - 10} more\n")
            
            # Tissue/Cell type
            f.write("\n\nTISSUE/CELL TYPE\n")
            f.write("-"*80 + "\n")
            tissue_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                          ['tissue', 'body_site', 'cell_type', 'cell_line', 'source_name'])]
            
            for col in tissue_cols[:5]:
                if col in self.df.columns:
                    non_null = self.df[col].dropna()
                    if len(non_null) > 0:
                        f.write(f"\n{col}:\n")
                        counts = self.df[col].value_counts().head(10)
                        for val, count in counts.items():
                            val_str = str(val)[:60]
                            f.write(f"  {val_str}: {count}\n")
            
            # Sequencing information
            f.write("\n\nSEQUENCING INFORMATION\n")
            f.write("-"*80 + "\n")
            
            if 'platform_type' in self.df.columns:
                f.write("\nPlatform:\n")
                counts = self.df['platform_type'].value_counts()
                for val, count in counts.items():
                    f.write(f"  {val}: {count}\n")
            
            if 'instrument_model' in self.df.columns:
                f.write("\nInstrument:\n")
                counts = self.df['instrument_model'].value_counts().head(10)
                for val, count in counts.items():
                    f.write(f"  {val}: {count}\n")
            
            if 'library_strategy' in self.df.columns:
                f.write("\nLibrary Strategy:\n")
                counts = self.df['library_strategy'].value_counts()
                for val, count in counts.items():
                    f.write(f"  {val}: {count}\n")
            
            if 'library_layout' in self.df.columns:
                f.write("\nLibrary Layout:\n")
                counts = self.df['library_layout'].value_counts()
                for val, count in counts.items():
                    f.write(f"  {val}: {count}\n")
            
            # Sequencing depth
            if 'run_total_spots' in self.df.columns:
                f.write("\nSequencing Depth (total spots):\n")
                spots = pd.to_numeric(self.df['run_total_spots'], errors='coerce').dropna()
                if len(spots) > 0:
                    f.write(f"  Mean: {spots.mean():,.0f}\n")
                    f.write(f"  Median: {spots.median():,.0f}\n")
                    f.write(f"  Min: {spots.min():,.0f}\n")
                    f.write(f"  Max: {spots.max():,.0f}\n")
            
            # Data completeness
            f.write("\n\nDATA COMPLETENESS\n")
            f.write("-"*80 + "\n")
            
            key_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                       ['sex', 'age', 'disease', 'tissue', 'body_site', 'phenotype', 
                        'treatment', 'genotype'])]
            
            if key_cols:
                completeness = []
                for col in key_cols:
                    pct = (self.df[col].notna().sum() / len(self.df)) * 100
                    completeness.append((col, pct, self.df[col].notna().sum()))
                
                completeness.sort(key=lambda x: x[1], reverse=True)
                
                for col, pct, count in completeness:
                    f.write(f"{col:50s} {pct:6.1f}% ({count}/{len(self.df)})\n")
            
            # Study information
            f.write("\n\nSTUDY INFORMATION\n")
            f.write("-"*80 + "\n")
            
            if 'study_accession' in self.df.columns:
                n_studies = self.df['study_accession'].nunique()
                f.write(f"Total unique studies: {n_studies}\n\n")
                f.write("Top studies by sample count:\n")
                counts = self.df['study_accession'].value_counts().head(10)
                for study, count in counts.items():
                    f.write(f"  {study}: {count} samples\n")
        
        print(f"Report saved to: {report_file}")
        return report_file
    
    def plot_demographics(self):
        """Plot demographic distributions."""
        
        # Sex distribution
        sex_cols = [c for c in self.sample_cols if any(x in c.lower() for x in ['sex', 'gender'])]
        if sex_cols:
            col = sex_cols[0]
            if self.df[col].notna().sum() > 0:
                plt.figure(figsize=(10, 6))
                counts = self.df[col].value_counts()
                plt.bar(counts.index, counts.values)
                plt.title('Sex/Gender Distribution')
                plt.xlabel('Sex/Gender')
                plt.ylabel('Count')
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(self.output_dir / 'demographics_sex.png', dpi=300)
                plt.close()
                print(f"Saved: demographics_sex.png")
    
    def plot_disease_distribution(self):
        """Plot disease distribution."""
        
        disease_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                       ['disease', 'phenotype', 'condition'])]
        
        for col in disease_cols[:2]:
            if col in self.df.columns and self.df[col].notna().sum() > 0:
                plt.figure(figsize=(12, 8))
                counts = self.df[col].value_counts().head(15)
                
                plt.barh(range(len(counts)), counts.values)
                plt.yticks(range(len(counts)), counts.index)
                plt.xlabel('Count')
                plt.title(f'Distribution: {col.replace("sample_", "")}')
                plt.tight_layout()
                
                filename = f'disease_{col.replace("sample_", "")}.png'
                plt.savefig(self.output_dir / filename, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Saved: {filename}")
    
    def plot_tissue_distribution(self):
        """Plot tissue/cell type distribution."""
        
        tissue_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                      ['tissue', 'body_site', 'cell_type'])]
        
        for col in tissue_cols[:2]:
            if col in self.df.columns and self.df[col].notna().sum() > 0:
                plt.figure(figsize=(12, 8))
                counts = self.df[col].value_counts().head(15)
                
                plt.barh(range(len(counts)), counts.values)
                plt.yticks(range(len(counts)), counts.index)
                plt.xlabel('Count')
                plt.title(f'Distribution: {col.replace("sample_", "")}')
                plt.tight_layout()
                
                filename = f'tissue_{col.replace("sample_", "")}.png'
                plt.savefig(self.output_dir / filename, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"Saved: {filename}")
    
    def plot_platform_distribution(self):
        """Plot sequencing platform distribution."""
        
        if 'platform_type' in self.df.columns:
            plt.figure(figsize=(10, 6))
            counts = self.df['platform_type'].value_counts()
            plt.pie(counts.values, labels=counts.index, autopct='%1.1f%%')
            plt.title('Sequencing Platform Distribution')
            plt.tight_layout()
            plt.savefig(self.output_dir / 'platform_distribution.png', dpi=300)
            plt.close()
            print(f"Saved: platform_distribution.png")
        
        if 'library_strategy' in self.df.columns:
            plt.figure(figsize=(10, 6))
            counts = self.df['library_strategy'].value_counts()
            plt.bar(counts.index, counts.values)
            plt.title('Library Strategy Distribution')
            plt.xlabel('Strategy')
            plt.ylabel('Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(self.output_dir / 'library_strategy.png', dpi=300)
            plt.close()
            print(f"Saved: library_strategy.png")
    
    def plot_completeness(self):
        """Plot data completeness."""
        
        key_cols = [c for c in self.sample_cols if any(x in c.lower() for x in 
                   ['sex', 'age', 'disease', 'tissue', 'body_site', 'phenotype'])][:20]
        
        if key_cols:
            completeness = []
            labels = []
            for col in key_cols:
                pct = (self.df[col].notna().sum() / len(self.df)) * 100
                completeness.append(pct)
                labels.append(col.replace('sample_', '')[:30])
            
            plt.figure(figsize=(12, 8))
            plt.barh(range(len(completeness)), completeness)
            plt.yticks(range(len(completeness)), labels)
            plt.xlabel('Completeness (%)')
            plt.title('Data Completeness for Key Patient Attributes')
            plt.xlim(0, 100)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            plt.savefig(self.output_dir / 'data_completeness.png', dpi=300)
            plt.close()
            print(f"Saved: data_completeness.png")
    
    def plot_sequencing_depth(self):
        """Plot sequencing depth distribution."""
        
        if 'run_total_spots' in self.df.columns:
            spots = pd.to_numeric(self.df['run_total_spots'], errors='coerce').dropna()
            
            if len(spots) > 0:
                plt.figure(figsize=(12, 6))
                
                plt.subplot(1, 2, 1)
                plt.hist(spots / 1e6, bins=50, edgecolor='black')
                plt.xlabel('Total Spots (millions)')
                plt.ylabel('Frequency')
                plt.title('Sequencing Depth Distribution')
                
                plt.subplot(1, 2, 2)
                plt.boxplot(spots / 1e6)
                plt.ylabel('Total Spots (millions)')
                plt.title('Sequencing Depth (Box Plot)')
                
                plt.tight_layout()
                plt.savefig(self.output_dir / 'sequencing_depth.png', dpi=300)
                plt.close()
                print(f"Saved: sequencing_depth.png")
    
    def run_full_analysis(self):
        """Run complete EDA pipeline."""
        
        print("\n" + "="*70)
        print("RUNNING METADATA EDA ANALYSIS")
        print("="*70)
        print(f"Output directory: {self.output_dir}")
        print()
        
        # Generate text report
        print("Generating text report...")
        report_file = self.generate_report()
        
        # Generate plots
        print("\nGenerating plots...")
        self.plot_demographics()
        self.plot_disease_distribution()
        self.plot_tissue_distribution()
        self.plot_platform_distribution()
        self.plot_completeness()
        self.plot_sequencing_depth()
        
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE")
        print("="*70)
        print(f"Report: {report_file}")
        print(f"Plots: {self.output_dir}/")
        print("\nGenerated files:")
        for f in sorted(self.output_dir.glob('*')):
            print(f"  - {f.name}")
        print("="*70)


def main():
    """Main entry point."""
    
    if len(sys.argv) < 2:
        # Find most recent metadata file
        metadata_dir = Path(__file__).parent.parent / "data" / "metadata"
        csv_files = list(metadata_dir.glob("sra_metadata_complete_*.csv"))
        
        if not csv_files:
            print("ERROR: No metadata files found!")
            print("Run: python scripts/download_all_metadata.py first")
            sys.exit(1)
        
        csv_file = max(csv_files, key=lambda p: p.stat().st_mtime)
        print(f"Using most recent file: {csv_file.name}\n")
    else:
        csv_file = sys.argv[1]
    
    # Run EDA
    eda = MetadataEDA(csv_file)
    eda.run_full_analysis()


if __name__ == "__main__":
    main()
