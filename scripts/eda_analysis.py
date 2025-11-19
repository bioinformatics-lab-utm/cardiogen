"""
Exploratory Data Analysis (EDA) for Cardiogen Pipeline Results
===============================================================

Analyzes QC results and generates visualizations and statistics.

Usage:
    python scripts/eda_analysis.py
    python scripts/eda_analysis.py --output-dir reports/
"""

import sys
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List

# Add to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config.config import RESULTS_DIR, DATA_DIR

# Try to import visualization libraries
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("⚠️  Warning: matplotlib/seaborn not installed. Plots will be skipped.")
    print("   Install with: pip install matplotlib seaborn")


class CardioEDA:
    """Exploratory Data Analysis for Cardiogen results"""
    
    def __init__(self, results_dir: Path, output_dir: Path):
        self.results_dir = results_dir
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set plotting style
        if PLOTTING_AVAILABLE:
            sns.set_style("whitegrid")
            plt.rcParams['figure.figsize'] = (12, 8)
    
    def load_data(self) -> pd.DataFrame:
        """Load processing summary data"""
        summary_file = self.results_dir / 'processing_summary.csv'
        
        if not summary_file.exists():
            raise FileNotFoundError(
                f"Summary file not found: {summary_file}\n"
                "Run the pipeline first: python scripts/run_pipeline.py"
            )
        
        df = pd.read_csv(summary_file)
        print(f"✅ Loaded {len(df)} samples from {summary_file}")
        return df
    
    def load_batch_results(self) -> List[Dict]:
        """Load all batch results"""
        batch_files = sorted(self.results_dir.glob('batch_*_results.json'))
        all_results = []
        
        for batch_file in batch_files:
            with open(batch_file) as f:
                all_results.extend(json.load(f))
        
        print(f"✅ Loaded {len(all_results)} results from {len(batch_files)} batches")
        return all_results
    
    def basic_statistics(self, df: pd.DataFrame) -> Dict:
        """Calculate basic statistics"""
        stats = {
            'total_samples': len(df),
            'diseases': df['disease'].nunique(),
            'total_sequences': df['total_sequences'].sum(),
            'avg_sequences_per_sample': df['total_sequences'].mean(),
            'total_data_mb': df['fastq_size_mb'].sum(),
            'avg_gc_content': df['gc_content'].mean(),
            'gc_std': df['gc_content'].std(),
        }
        
        return stats
    
    def disease_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Statistics grouped by disease"""
        disease_stats = df.groupby('disease').agg({
            'accession': 'count',
            'total_sequences': ['sum', 'mean', 'std'],
            'gc_content': ['mean', 'std'],
            'fastq_size_mb': 'sum',
            'sequence_length': lambda x: x.mode()[0] if len(x) > 0 else None
        }).round(2)
        
        disease_stats.columns = [
            'n_samples', 
            'total_sequences', 'avg_sequences', 'std_sequences',
            'avg_gc', 'std_gc',
            'total_data_mb',
            'common_length'
        ]
        
        return disease_stats
    
    def generate_report(self, df: pd.DataFrame, stats: Dict, 
                       disease_stats: pd.DataFrame) -> None:
        """Generate text report"""
        report_file = self.output_dir / 'eda_report.txt'
        
        report = f"""
{'='*70}
Cardiogen EDA Report
{'='*70}
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

OVERALL STATISTICS
{'='*70}
Total Samples:              {stats['total_samples']:,}
Unique Diseases:            {stats['diseases']}
Total Sequences:            {stats['total_sequences']:,.0f}
Avg Sequences/Sample:       {stats['avg_sequences_per_sample']:,.0f}
Total Data Processed:       {stats['total_data_mb']:,.2f} MB ({stats['total_data_mb']/1024:.2f} GB)
Average GC Content:         {stats['avg_gc_content']:.2f}% (± {stats['gc_std']:.2f})

STATISTICS BY DISEASE
{'='*70}
{disease_stats.to_string()}

DATA QUALITY INDICATORS
{'='*70}
GC Content Range:           {df['gc_content'].min():.1f}% - {df['gc_content'].max():.1f}%
Sequence Length Range:      {df['sequence_length'].min()} - {df['sequence_length'].max()}
Most Common Length:         {df['sequence_length'].mode()[0] if len(df) > 0 else 'N/A'}

Samples by Disease:
{df['disease'].value_counts().to_string()}

POTENTIAL QUALITY ISSUES
{'='*70}
"""
        
        # Check for outliers
        gc_mean = df['gc_content'].mean()
        gc_std = df['gc_content'].std()
        outliers_gc = df[(df['gc_content'] < gc_mean - 2*gc_std) | 
                         (df['gc_content'] > gc_mean + 2*gc_std)]
        
        if len(outliers_gc) > 0:
            report += f"⚠️  {len(outliers_gc)} samples with unusual GC content (>2 SD from mean)\n"
            report += f"   Accessions: {', '.join(outliers_gc['accession'].head(10).tolist())}\n"
        else:
            report += "✅ No significant GC content outliers detected\n"
        
        # Check sequence counts
        seq_mean = df['total_sequences'].mean()
        low_seq = df[df['total_sequences'] < seq_mean * 0.5]
        
        if len(low_seq) > 0:
            report += f"\n⚠️  {len(low_seq)} samples with low sequence count (<50% of mean)\n"
            report += f"   Accessions: {', '.join(low_seq['accession'].head(10).tolist())}\n"
        else:
            report += "\n✅ All samples have adequate sequence counts\n"
        
        report += f"\n{'='*70}\n"
        report += "End of Report\n"
        report += f"{'='*70}\n"
        
        with open(report_file, 'w') as f:
            f.write(report)
        
        print(f"✅ Report saved to: {report_file}")
        print(report)
    
    def plot_gc_distribution(self, df: pd.DataFrame) -> None:
        """Plot GC content distribution"""
        if not PLOTTING_AVAILABLE:
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Overall distribution
        axes[0].hist(df['gc_content'], bins=30, edgecolor='black', alpha=0.7)
        axes[0].axvline(df['gc_content'].mean(), color='red', 
                       linestyle='--', label=f'Mean: {df["gc_content"].mean():.1f}%')
        axes[0].set_xlabel('GC Content (%)')
        axes[0].set_ylabel('Frequency')
        axes[0].set_title('GC Content Distribution - All Samples')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # By disease
        diseases = df['disease'].unique()
        for disease in diseases:
            disease_data = df[df['disease'] == disease]['gc_content']
            axes[1].hist(disease_data, bins=20, alpha=0.5, label=disease)
        
        axes[1].set_xlabel('GC Content (%)')
        axes[1].set_ylabel('Frequency')
        axes[1].set_title('GC Content Distribution by Disease')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = self.output_dir / 'gc_content_distribution.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved: {output_file}")
    
    def plot_gc_boxplot(self, df: pd.DataFrame) -> None:
        """Boxplot of GC content by disease"""
        if not PLOTTING_AVAILABLE:
            return
        
        plt.figure(figsize=(12, 6))
        df_sorted = df.sort_values('disease')
        sns.boxplot(data=df_sorted, x='disease', y='gc_content', palette='Set2')
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Disease')
        plt.ylabel('GC Content (%)')
        plt.title('GC Content Distribution by Disease')
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        
        output_file = self.output_dir / 'gc_content_boxplot.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved: {output_file}")
    
    def plot_sequence_counts(self, df: pd.DataFrame) -> None:
        """Plot sequence counts"""
        if not PLOTTING_AVAILABLE:
            return
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # Bar plot by disease
        disease_counts = df.groupby('disease')['total_sequences'].sum() / 1e6
        disease_counts.plot(kind='bar', ax=axes[0], color='steelblue')
        axes[0].set_xlabel('Disease')
        axes[0].set_ylabel('Total Sequences (Millions)')
        axes[0].set_title('Total Sequences by Disease')
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].grid(True, alpha=0.3, axis='y')
        
        # Violin plot
        sns.violinplot(data=df, x='disease', y='total_sequences', ax=axes[1])
        axes[1].set_xlabel('Disease')
        axes[1].set_ylabel('Sequences per Sample')
        axes[1].set_title('Sequence Count Distribution by Disease')
        axes[1].tick_params(axis='x', rotation=45)
        axes[1].set_yscale('log')
        
        plt.tight_layout()
        output_file = self.output_dir / 'sequence_counts.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved: {output_file}")
    
    def plot_samples_per_disease(self, df: pd.DataFrame) -> None:
        """Plot number of samples per disease"""
        if not PLOTTING_AVAILABLE:
            return
        
        plt.figure(figsize=(10, 6))
        disease_counts = df['disease'].value_counts()
        disease_counts.plot(kind='barh', color='coral')
        plt.xlabel('Number of Samples')
        plt.ylabel('Disease')
        plt.title('Sample Distribution by Disease')
        plt.grid(True, alpha=0.3, axis='x')
        plt.tight_layout()
        
        output_file = self.output_dir / 'samples_per_disease.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved: {output_file}")
    
    def plot_correlation_matrix(self, df: pd.DataFrame) -> None:
        """Plot correlation matrix of numeric features"""
        if not PLOTTING_AVAILABLE:
            return
        
        # Select numeric columns
        numeric_cols = ['total_sequences', 'gc_content', 'fastq_size_mb']
        numeric_df = df[numeric_cols].dropna()
        
        if len(numeric_df) == 0:
            return
        
        plt.figure(figsize=(8, 6))
        correlation = numeric_df.corr()
        sns.heatmap(correlation, annot=True, cmap='coolwarm', center=0,
                   square=True, linewidths=1, cbar_kws={"shrink": 0.8})
        plt.title('Correlation Matrix of Numeric Features')
        plt.tight_layout()
        
        output_file = self.output_dir / 'correlation_matrix.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Plot saved: {output_file}")
    
    def run_full_analysis(self) -> None:
        """Run complete EDA analysis"""
        print("\n" + "="*70)
        print("🔬 Starting Cardiogen EDA Analysis")
        print("="*70 + "\n")
        
        # Load data
        df = self.load_data()
        
        # Calculate statistics
        print("\n📊 Calculating statistics...")
        stats = self.basic_statistics(df)
        disease_stats = self.disease_statistics(df)
        
        # Generate report
        print("\n📝 Generating report...")
        self.generate_report(df, stats, disease_stats)
        
        # Save detailed CSV
        csv_file = self.output_dir / 'disease_statistics.csv'
        disease_stats.to_csv(csv_file)
        print(f"✅ Detailed stats saved: {csv_file}")
        
        # Generate plots
        if PLOTTING_AVAILABLE:
            print("\n📈 Generating visualizations...")
            self.plot_samples_per_disease(df)
            self.plot_gc_distribution(df)
            self.plot_gc_boxplot(df)
            self.plot_sequence_counts(df)
            self.plot_correlation_matrix(df)
        
        print("\n" + "="*70)
        print("✅ EDA Analysis Complete!")
        print(f"📁 Results saved in: {self.output_dir}")
        print("="*70 + "\n")
        
        # List output files
        print("Generated files:")
        for file in sorted(self.output_dir.glob('*')):
            print(f"  - {file.name}")


def main(output_dir: str = None):
    """Main function"""
    if output_dir is None:
        output_dir = RESULTS_DIR / 'eda'
    else:
        output_dir = Path(output_dir)
    
    eda = CardioEDA(
        results_dir=RESULTS_DIR,
        output_dir=output_dir
    )
    
    try:
        eda.run_full_analysis()
    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        print("\nPlease run the pipeline first:")
        print("  python scripts/run_pipeline.py --batch-size 5")
        return 1
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Exploratory Data Analysis for Cardiogen results'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory for EDA results (default: data/results/eda/)'
    )
    
    args = parser.parse_args()
    exit(main(args.output_dir))
