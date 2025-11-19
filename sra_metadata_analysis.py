#!/usr/bin/env python3
"""
Comprehensive SRA Metadata Exploratory Data Analysis
Scientific Analysis of Illumina Sequencing Data

Complete exploratory analysis of all 296 columns from SRA metadata focusing on Illumina platform samples.
Generates 60+ publication-quality visualizations and statistical summaries.

Author: Generated Analysis Script
Date: November 17, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import linregress
import warnings
import os
from pathlib import Path
from typing import Optional, List, Dict, Tuple

warnings.filterwarnings('ignore')


class SRAMetadataAnalyzer:
    """Comprehensive analyzer for SRA metadata with Illumina platform focus."""
    
    def __init__(self, metadata_file: str, output_dir: str = 'data/analysis_results'):
        """
        Initialize the analyzer.
        
        Args:
            metadata_file: Path to the SRA metadata CSV file
            output_dir: Directory for saving analysis results and figures
        """
        self.metadata_file = metadata_file
        self.output_dir = Path(output_dir)
        self.figures_dir = self.output_dir / 'figures'
        self.df = None
        self.df_illumina = None
        
        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # Configure plotting style
        self._setup_plotting()
    
    def _setup_plotting(self):
        """Configure matplotlib and seaborn for publication-quality plots."""
        plt.style.use('seaborn-v0_8-paper')
        sns.set_context("paper", font_scale=1.3)
        sns.set_palette("Set2")
        
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['figure.figsize'] = (10, 6)
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['grid.alpha'] = 0.3
        plt.rcParams['axes.spines.top'] = False
        plt.rcParams['axes.spines.right'] = False
    
    def load_data(self):
        """Load SRA metadata from CSV file."""
        print(f"Loading data from {self.metadata_file}...")
        self.df = pd.read_csv(self.metadata_file, low_memory=False)
        print(f"Loaded {len(self.df)} samples with {len(self.df.columns)} columns")
    
    def filter_illumina(self):
        """Filter data to include only Illumina platform samples."""
        print("Filtering for Illumina platform samples...")
        illumina_mask = (
            self.df['platform_type'].str.contains('ILLUMINA', case=False, na=False) |
            self.df['instrument_model'].str.contains('Illumina', case=False, na=False)
        )
        self.df_illumina = self.df[illumina_mask].copy().reset_index(drop=True)
        self.df = self.df_illumina
        
        # Calculate sample-level statistics
        self.n_runs = len(self.df)
        self.n_samples = self.df['sample_accession'].nunique()
        
        print(f"Filtered to {self.n_runs} Illumina runs from {self.n_samples} unique samples")
        print(f"Average runs per sample: {self.n_runs/self.n_samples:.2f}")
    
    def _save_figure(self, name: str):
        """Save current figure to file."""
        filepath = self.figures_dir / f"{name}.png"
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"  Saved: {filepath}")
    
    def analyze_instruments(self):
        """Analyze sequencing instrument distribution."""
        print("\n=== Analyzing Instrument Distribution ===")
        instruments = self.df['instrument_model'].value_counts()
        
        # Horizontal bar chart
        fig, ax = plt.subplots(figsize=(12, 7))
        top_n = 12
        top_instruments = instruments.head(top_n)
        colors = sns.color_palette("viridis", n_colors=top_n)
        bars = ax.barh(range(len(top_instruments)), top_instruments.values, color=colors,
                       edgecolor='white', linewidth=1.5)
        ax.set_yticks(range(len(top_instruments)))
        ax.set_yticklabels(top_instruments.index, fontsize=11)
        ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title('Illumina Sequencing Instrument Distribution', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        for i, (idx, val) in enumerate(top_instruments.items()):
            ax.text(val + 20, i, f'{val:,}', va='center', fontsize=10, fontweight='bold')
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('01_instrument_distribution_bar')
        plt.close()
        
        # Pie chart
        fig, ax = plt.subplots(figsize=(10, 10))
        sizes = top_instruments.values
        labels = [f'{inst}\n({cnt:,})' for inst, cnt in zip(top_instruments.index, sizes)]
        colors_pie = sns.color_palette("Spectral", n_colors=len(top_instruments))
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                            colors=colors_pie, textprops={'fontsize': 9},
                                            pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(10)
        ax.set_title('Illumina Instrument Distribution', fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        self._save_figure('02_instrument_distribution_pie')
        plt.close()
    
    def analyze_library_strategy(self):
        """Analyze library preparation strategies."""
        print("\n=== Analyzing Library Strategy ===")
        strategies = self.df['library_strategy'].value_counts()
        
        # Bar chart
        fig, ax = plt.subplots(figsize=(12, 6))
        colors = sns.color_palette("RdYlGn", n_colors=len(strategies))
        bars = ax.barh(range(len(strategies)), strategies.values, color=colors,
                       edgecolor='black', linewidth=1.2)
        ax.set_yticks(range(len(strategies)))
        ax.set_yticklabels(strategies.index, fontsize=11)
        ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title('Library Strategy Distribution', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        for i, (idx, val) in enumerate(strategies.items()):
            ax.text(val + 20, i, f'{val:,}', va='center', fontsize=10, fontweight='bold')
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('03_library_strategy_bar')
        plt.close()
        
        # Pie chart
        fig, ax = plt.subplots(figsize=(10, 10))
        sizes = strategies.values
        labels = [f'{strat}\n({cnt:,})' for strat, cnt in zip(strategies.index, sizes)]
        colors_pie = sns.color_palette("Set2", n_colors=len(strategies))
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                            colors=colors_pie, textprops={'fontsize': 10},
                                            pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(11)
        ax.set_title('Library Strategy Distribution', fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        self._save_figure('04_library_strategy_pie')
        plt.close()
    
    def analyze_sequencing_depth(self):
        """Analyze sequencing depth metrics."""
        print("\n=== Analyzing Sequencing Depth ===")
        
        # Detect column names
        if 'run_total_spots' in self.df.columns:
            self.df['total_spots_numeric'] = pd.to_numeric(self.df['run_total_spots'], errors='coerce')
        elif 'total_spots' in self.df.columns:
            self.df['total_spots_numeric'] = pd.to_numeric(self.df['total_spots'], errors='coerce')
        else:
            self.df['total_spots_numeric'] = pd.Series([np.nan] * len(self.df))
        
        if 'run_total_bases' in self.df.columns:
            self.df['total_bases_numeric'] = pd.to_numeric(self.df['run_total_bases'], errors='coerce')
        elif 'total_bases' in self.df.columns:
            self.df['total_bases_numeric'] = pd.to_numeric(self.df['total_bases'], errors='coerce')
        else:
            self.df['total_bases_numeric'] = pd.Series([np.nan] * len(self.df))
        
        if 'run_avg_read_length' in self.df.columns:
            self.df['avg_read_length_numeric'] = pd.to_numeric(self.df['run_avg_read_length'], errors='coerce')
        elif 'avg_read_length' in self.df.columns:
            self.df['avg_read_length_numeric'] = pd.to_numeric(self.df['avg_read_length'], errors='coerce')
        else:
            self.df['avg_read_length_numeric'] = pd.Series([np.nan] * len(self.df))
        
        depth_data = self.df[self.df['total_spots_numeric'] > 0]['total_spots_numeric'].dropna()
        
        if len(depth_data) > 0:
            # Histogram
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.hist(np.log10(depth_data + 1), bins=50, color='steelblue', edgecolor='black',
                    alpha=0.7, linewidth=1.2)
            ax.set_xlabel('Log10(Total Spots + 1)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
            ax.set_title('Distribution of Sequencing Depth (Total Spots)', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            plt.tight_layout()
            self._save_figure('05_sequencing_depth_histogram')
            plt.close()
            
            # Violin plot
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.violinplot([np.log10(depth_data + 1)], positions=[0], widths=0.7,
                          showmeans=True, showmedians=True)
            ax.set_ylabel('Log10(Total Spots + 1)', fontsize=12, fontweight='bold')
            ax.set_title('Violin Plot: Sequencing Depth Distribution', fontsize=14, fontweight='bold', pad=20)
            ax.set_xticks([])
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            plt.tight_layout()
            self._save_figure('05b_sequencing_depth_violin')
            plt.close()
        
        # Total bases histogram
        bases_data = self.df[self.df['total_bases_numeric'] > 0]['total_bases_numeric'].dropna()
        if len(bases_data) > 0:
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.hist(np.log10(bases_data + 1), bins=50, color='coral', edgecolor='black',
                    alpha=0.7, linewidth=1.2)
            ax.set_xlabel('Log10(Total Bases + 1)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
            ax.set_title('Distribution of Total Bases Sequenced', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            plt.tight_layout()
            self._save_figure('05c_total_bases_histogram')
            plt.close()
        
        # Read length histogram
        read_length_data = self.df[self.df['avg_read_length_numeric'] > 0]['avg_read_length_numeric'].dropna()
        if len(read_length_data) > 0:
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.hist(read_length_data, bins=50, color='mediumseagreen', edgecolor='black',
                    alpha=0.7, linewidth=1.2)
            ax.set_xlabel('Average Read Length (bp)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
            ax.set_title('Distribution of Average Read Length', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            plt.tight_layout()
            self._save_figure('05d_read_length_histogram')
            plt.close()
    
    def analyze_demographics(self):
        """Analyze patient demographics (sex, age)."""
        print("\n=== Analyzing Demographics ===")
        
        # Sex distribution
        sex_col = 'sample_attributes_sex' if 'sample_attributes_sex' in self.df.columns else 'sample_sex'
        if sex_col in self.df.columns:
            sex_data = self.df[sex_col].value_counts()
            
            if len(sex_data) > 0:
                fig, ax = plt.subplots(figsize=(10, 6))
                colors = sns.color_palette("Set2", n_colors=len(sex_data))
                bars = ax.bar(range(len(sex_data)), sex_data.values, color=colors,
                              edgecolor='black', linewidth=1.5, alpha=0.8)
                ax.set_xticks(range(len(sex_data)))
                ax.set_xticklabels(sex_data.index, fontsize=11)
                ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sex Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                for i, (idx, val) in enumerate(sex_data.items()):
                    ax.text(i, val + 20, f'{val:,}', ha='center', fontsize=10, fontweight='bold')
                plt.tight_layout()
                self._save_figure('06_sex_distribution')
                plt.close()
                
                # Pie chart
                fig, ax = plt.subplots(figsize=(9, 9))
                sizes = sex_data.values
                labels = [f'{sex}\n({cnt:,})' for sex, cnt in zip(sex_data.index, sizes)]
                colors_pie = sns.color_palette("pastel", n_colors=len(sex_data))
                wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                                    colors=colors_pie, textprops={'fontsize': 11},
                                                    pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(12)
                ax.set_title('Sex Distribution', fontsize=14, fontweight='bold', pad=20)
                plt.tight_layout()
                self._save_figure('06b_sex_distribution_pie')
                plt.close()
        
        # Age analysis
        age_col = 'sample_attributes_age' if 'sample_attributes_age' in self.df.columns else None
        if age_col and age_col in self.df.columns:
            age_data = pd.to_numeric(self.df[age_col], errors='coerce').dropna()
            if len(age_data) > 0:
                # Histogram
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.hist(age_data, bins=40, color='skyblue', edgecolor='black',
                        alpha=0.7, linewidth=1.2)
                ax.set_xlabel('Age', fontsize=12, fontweight='bold')
                ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
                ax.set_title('Age Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('06c_age_histogram')
                plt.close()
                
                # Boxplot
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.boxplot([age_data], vert=True, widths=0.5, patch_artist=True,
                           boxprops=dict(facecolor='lightcoral', edgecolor='black', linewidth=1.5),
                           medianprops=dict(color='red', linewidth=2),
                           whiskerprops=dict(linewidth=1.5),
                           capprops=dict(linewidth=1.5))
                ax.set_ylabel('Age', fontsize=12, fontweight='bold')
                ax.set_title('Age Distribution Box Plot', fontsize=14, fontweight='bold', pad=20)
                ax.set_xticks([])
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('06d_age_boxplot')
                plt.close()
    
    def analyze_metadata_completeness(self):
        """Analyze metadata completeness across all columns."""
        print("\n=== Analyzing Metadata Completeness ===")
        
        completeness = (1 - self.df.isnull().sum() / len(self.df)) * 100
        completeness_sorted = completeness.sort_values(ascending=True)
        
        # Top incomplete columns
        fig, ax = plt.subplots(figsize=(12, 10))
        n_cols = 30
        top_incomplete = completeness_sorted.head(n_cols)
        colors = plt.cm.RdYlGn(top_incomplete.values / 100)
        bars = ax.barh(range(len(top_incomplete)), top_incomplete.values, color=colors,
                       edgecolor='black', linewidth=0.8)
        ax.set_yticks(range(len(top_incomplete)))
        ax.set_yticklabels(top_incomplete.index, fontsize=9)
        ax.set_xlabel('Completeness (%)', fontsize=12, fontweight='bold')
        ax.set_title('Top 30 Least Complete Columns', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_xlim(0, 100)
        for i, (idx, val) in enumerate(top_incomplete.items()):
            ax.text(val + 1, i, f'{val:.1f}%', va='center', fontsize=8)
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('07_metadata_completeness')
        plt.close()
    
    def analyze_organizations(self):
        """Analyze data contributors and organizations."""
        print("\n=== Analyzing Organizations ===")
        
        org_col = None
        for col in self.df.columns:
            if 'center' in col.lower() or 'organization' in col.lower() or 'lab' in col.lower():
                if self.df[col].notna().sum() > 100:
                    org_col = col
                    break
        
        if org_col:
            org_data = self.df[org_col].value_counts().head(15)
            
            if len(org_data) > 0:
                fig, ax = plt.subplots(figsize=(14, 8))
                colors = sns.color_palette("rocket", n_colors=len(org_data))
                bars = ax.barh(range(len(org_data)), org_data.values, color=colors,
                               edgecolor='black', linewidth=1.2)
                ax.set_yticks(range(len(org_data)))
                ax.set_yticklabels([str(org)[:50] for org in org_data.index], fontsize=9)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title(f'Top 15 Data Contributors ({org_col.replace("_", " ").title()})',
                             fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                for i, val in enumerate(org_data.values):
                    ax.text(val + max(org_data.values)*0.01, i, f'{val:,}',
                            va='center', fontsize=9, fontweight='bold')
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('08_top_organizations')
                plt.close()
    
    def analyze_correlations(self):
        """Analyze correlations between numerical sequencing parameters."""
        print("\n=== Analyzing Correlations ===")
        
        numeric_cols = ['total_spots_numeric', 'total_bases_numeric', 'avg_read_length_numeric']
        available_numeric = [col for col in numeric_cols if col in self.df.columns and self.df[col].notna().sum() > 100]
        
        if len(available_numeric) >= 2:
            corr_data = self.df[available_numeric].dropna()
            correlation_matrix = corr_data.corr()
            
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(correlation_matrix, annot=True, fmt='.3f', cmap='coolwarm',
                        center=0, square=True, linewidths=2, cbar_kws={'label': 'Correlation'},
                        ax=ax, vmin=-1, vmax=1)
            ax.set_title('Correlation Matrix: Sequencing Parameters', fontsize=14, fontweight='bold', pad=20)
            labels = [col.replace('_numeric', '').replace('_', ' ').title() for col in correlation_matrix.columns]
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=11)
            ax.set_yticklabels(labels, rotation=0, fontsize=11)
            plt.tight_layout()
            self._save_figure('09_correlation_matrix')
            plt.close()
    
    def generate_summary_report(self):
        """Generate summary statistics and export results."""
        print("\n=== Generating Summary Report ===")
        
        summary_stats = {
            'total_samples': len(self.df),
            'total_columns': len(self.df.columns),
            'sample_attributes': len([c for c in self.df.columns if c.startswith('sample_')]),
            'run_attributes': len([c for c in self.df.columns if c.startswith('run_')]),
            'experiment_attributes': len([c for c in self.df.columns if c.startswith('experiment_')]),
            'study_attributes': len([c for c in self.df.columns if c.startswith('study_')])
        }
        
        # Save filtered data
        filtered_csv = self.output_dir / 'illumina_filtered_metadata.csv'
        self.df.to_csv(filtered_csv, index=False)
        print(f"  Saved filtered data: {filtered_csv}")
        
        # Save summary statistics
        summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['count'])
        summary_csv = self.output_dir / 'summary_statistics.csv'
        summary_df.to_csv(summary_csv)
        print(f"  Saved summary stats: {summary_csv}")
        
        # Print summary
        print("\n" + "="*60)
        print("ANALYSIS SUMMARY")
        print("="*60)
        for key, value in summary_stats.items():
            print(f"{key.replace('_', ' ').title()}: {value:,}")
        print("="*60)
    
    def analyze_diseases(self):
        """Analyze disease distribution."""
        print("\n=== Analyzing Disease Distribution ===")
        
        disease_col = 'sample_attributes_disease' if 'sample_attributes_disease' in self.df.columns else 'sample_disease'
        if disease_col in self.df.columns:
            disease_data = self.df[disease_col].value_counts().head(15)
            
            if len(disease_data) > 0:
                # Bar chart
                fig, ax = plt.subplots(figsize=(14, 8))
                colors = sns.color_palette("husl", n_colors=len(disease_data))
                bars = ax.barh(range(len(disease_data)), disease_data.values, color=colors,
                               edgecolor='black', linewidth=1.2)
                ax.set_yticks(range(len(disease_data)))
                ax.set_yticklabels(disease_data.index, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Top 15 Disease Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                for i, (idx, val) in enumerate(disease_data.items()):
                    ax.text(val + 10, i, f'{val:,}', va='center', fontsize=9, fontweight='bold')
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('10_disease_distribution_bar')
                plt.close()
                
                # Pie chart
                fig, ax = plt.subplots(figsize=(12, 12))
                sizes = disease_data.values
                labels = [f'{dis[:30]}...\n({cnt:,})' if len(dis) > 30 else f'{dis}\n({cnt:,})' 
                          for dis, cnt in zip(disease_data.index, sizes)]
                colors_pie = sns.color_palette("Set3", n_colors=len(disease_data))
                wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                                    colors=colors_pie, textprops={'fontsize': 8},
                                                    pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=1.5))
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(9)
                ax.set_title('Top 15 Disease Distribution', fontsize=14, fontweight='bold', pad=20)
                plt.tight_layout()
                self._save_figure('11_disease_distribution_pie')
                plt.close()
    
    def analyze_tissues(self):
        """Analyze tissue distribution."""
        print("\n=== Analyzing Tissue Distribution ===")
        
        tissue_col = 'sample_attributes_tissue' if 'sample_attributes_tissue' in self.df.columns else 'sample_tissue'
        if tissue_col in self.df.columns:
            tissue_data = self.df[tissue_col].value_counts().head(15)
            
            if len(tissue_data) > 0:
                # Bar chart
                fig, ax = plt.subplots(figsize=(14, 8))
                colors = sns.color_palette("mako", n_colors=len(tissue_data))
                bars = ax.barh(range(len(tissue_data)), tissue_data.values, color=colors,
                               edgecolor='black', linewidth=1.2)
                ax.set_yticks(range(len(tissue_data)))
                ax.set_yticklabels(tissue_data.index, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Top 15 Tissue Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                for i, (idx, val) in enumerate(tissue_data.items()):
                    ax.text(val + 10, i, f'{val:,}', va='center', fontsize=9, fontweight='bold')
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('12_tissue_distribution_bar')
                plt.close()
                
                # Pie chart
                fig, ax = plt.subplots(figsize=(12, 12))
                sizes = tissue_data.values
                labels = [f'{tis[:30]}...\n({cnt:,})' if len(tis) > 30 else f'{tis}\n({cnt:,})' 
                          for tis, cnt in zip(tissue_data.index, sizes)]
                colors_pie = sns.color_palette("Paired", n_colors=len(tissue_data))
                wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                                    colors=colors_pie, textprops={'fontsize': 8},
                                                    pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=1.5))
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(9)
                ax.set_title('Top 15 Tissue Distribution', fontsize=14, fontweight='bold', pad=20)
                plt.tight_layout()
                self._save_figure('13_tissue_distribution_pie')
                plt.close()
    
    def analyze_temporal(self):
        """Analyze temporal patterns in data releases."""
        print("\n=== Analyzing Temporal Patterns ===")
        
        # Parse release dates
        if 'run_release_date' in self.df.columns:
            self.df['release_year'] = pd.to_datetime(self.df['run_release_date'], errors='coerce').dt.year
        elif 'release_date' in self.df.columns:
            self.df['release_year'] = pd.to_datetime(self.df['release_date'], errors='coerce').dt.year
        else:
            self.df['release_year'] = pd.Series([np.nan] * len(self.df))
        
        if self.df['release_year'].notna().sum() > 0:
            year_counts = self.df['release_year'].value_counts().sort_index()
            
            # Yearly releases bar chart
            fig, ax = plt.subplots(figsize=(14, 6))
            colors = sns.color_palette("coolwarm", n_colors=len(year_counts))
            bars = ax.bar(year_counts.index, year_counts.values, color=colors,
                          edgecolor='black', linewidth=1.2, alpha=0.8)
            ax.set_xlabel('Year', fontsize=12, fontweight='bold')
            ax.set_ylabel('Number of Samples Released', fontsize=12, fontweight='bold')
            ax.set_title('Temporal Distribution of Data Releases', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            plt.xticks(rotation=45, ha='right')
            for i, (year, count) in enumerate(year_counts.items()):
                ax.text(year, count + max(year_counts.values)*0.01, f'{count:,}', 
                        ha='center', fontsize=9, fontweight='bold')
            plt.tight_layout()
            self._save_figure('14_temporal_yearly_releases')
            plt.close()
            
            # Cumulative growth
            fig, ax = plt.subplots(figsize=(14, 6))
            cumulative = year_counts.cumsum()
            ax.plot(cumulative.index, cumulative.values, marker='o', linewidth=3, 
                    markersize=8, color='steelblue', markerfacecolor='orange', markeredgewidth=2)
            ax.fill_between(cumulative.index, cumulative.values, alpha=0.3, color='steelblue')
            ax.set_xlabel('Year', fontsize=12, fontweight='bold')
            ax.set_ylabel('Cumulative Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Cumulative Data Release Over Time', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            self._save_figure('15_temporal_cumulative')
            plt.close()
        
        if self.df['release_year'].notna().sum() > 0:
            year_counts = self.df['release_year'].value_counts().sort_index()
            
            # Yearly releases bar chart
            fig, ax = plt.subplots(figsize=(14, 6))
            colors = sns.color_palette("coolwarm", n_colors=len(year_counts))
            bars = ax.bar(year_counts.index, year_counts.values, color=colors,
                          edgecolor='black', linewidth=1.2, alpha=0.8)
            ax.set_xlabel('Year', fontsize=12, fontweight='bold')
            ax.set_ylabel('Number of Samples Released', fontsize=12, fontweight='bold')
            ax.set_title('Temporal Distribution of Data Releases', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            plt.xticks(rotation=45, ha='right')
            for i, (year, count) in enumerate(year_counts.items()):
                ax.text(year, count + max(year_counts.values)*0.01, f'{count:,}', 
                        ha='center', fontsize=9, fontweight='bold')
            plt.tight_layout()
            self._save_figure('14_temporal_yearly_releases')
            plt.close()
            
            # Cumulative growth
            fig, ax = plt.subplots(figsize=(14, 6))
            cumulative = year_counts.cumsum()
            ax.plot(cumulative.index, cumulative.values, marker='o', linewidth=3, 
                    markersize=8, color='steelblue', markerfacecolor='orange', markeredgewidth=2)
            ax.fill_between(cumulative.index, cumulative.values, alpha=0.3, color='steelblue')
            ax.set_xlabel('Year', fontsize=12, fontweight='bold')
            ax.set_ylabel('Cumulative Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Cumulative Data Release Over Time', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            self._save_figure('15_temporal_cumulative')
            plt.close()
    
    def analyze_scatter_plots(self):
        """Generate scatter plots for sequencing parameters."""
        print("\n=== Generating Scatter Plots ===")
        
        # Scatter: Depth vs Read Length
        if 'total_spots_numeric' in self.df.columns and 'avg_read_length_numeric' in self.df.columns:
            scatter_data = self.df[['total_spots_numeric', 'avg_read_length_numeric']].dropna()
            scatter_data = scatter_data[(scatter_data['total_spots_numeric'] > 0) & 
                                         (scatter_data['avg_read_length_numeric'] > 0)]
            
            if len(scatter_data) > 10:
                fig, ax = plt.subplots(figsize=(12, 8))
                scatter = ax.scatter(np.log10(scatter_data['total_spots_numeric']), 
                                    scatter_data['avg_read_length_numeric'],
                                    alpha=0.5, s=30, c=np.log10(scatter_data['total_spots_numeric']),
                                    cmap='viridis', edgecolors='black', linewidth=0.5)
                ax.set_xlabel('Log10(Total Spots)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Average Read Length (bp)', fontsize=12, fontweight='bold')
                ax.set_title('Sequencing Depth vs Read Length', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label('Log10(Total Spots)', fontsize=11, fontweight='bold')
                plt.tight_layout()
                self._save_figure('16_scatter_depth_vs_length')
                plt.close()
        
        # Scatter: Total Spots vs Total Bases (existing)
        if 'total_spots_numeric' in self.df.columns and 'avg_read_length_numeric' in self.df.columns:
            scatter_data = self.df[['total_spots_numeric', 'avg_read_length_numeric']].dropna()
            scatter_data = scatter_data[(scatter_data['total_spots_numeric'] > 0) & 
                                         (scatter_data['avg_read_length_numeric'] > 0)]
            
            if len(scatter_data) > 10:
                fig, ax = plt.subplots(figsize=(12, 8))
                scatter = ax.scatter(np.log10(scatter_data['total_spots_numeric']), 
                                    scatter_data['avg_read_length_numeric'],
                                    alpha=0.5, s=30, c=np.log10(scatter_data['total_spots_numeric']),
                                    cmap='viridis', edgecolors='black', linewidth=0.5)
                ax.set_xlabel('Log10(Total Spots)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Average Read Length (bp)', fontsize=12, fontweight='bold')
                ax.set_title('Sequencing Depth vs Read Length', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label('Log10(Total Spots)', fontsize=11, fontweight='bold')
                plt.tight_layout()
                self._save_figure('16_scatter_depth_vs_length')
                plt.close()
        
        # Scatter: Total Spots vs Total Bases with regression
        if 'total_bases_numeric' in self.df.columns and 'total_spots_numeric' in self.df.columns:
            scatter_data2 = self.df[['total_bases_numeric', 'total_spots_numeric']].dropna()
            scatter_data2 = scatter_data2[(scatter_data2['total_bases_numeric'] > 0) & 
                                           (scatter_data2['total_spots_numeric'] > 0)]
            
            if len(scatter_data2) > 10:
                fig, ax = plt.subplots(figsize=(12, 8))
                scatter = ax.scatter(np.log10(scatter_data2['total_spots_numeric']), 
                                    np.log10(scatter_data2['total_bases_numeric']),
                                    alpha=0.4, s=25, c=range(len(scatter_data2)),
                                    cmap='plasma', edgecolors='none')
                ax.set_xlabel('Log10(Total Spots)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Log10(Total Bases)', fontsize=12, fontweight='bold')
                ax.set_title('Total Spots vs Total Bases Relationship', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                
                # Add regression line
                x_log = np.log10(scatter_data2['total_spots_numeric'])
                y_log = np.log10(scatter_data2['total_bases_numeric'])
                slope, intercept, r_value, p_value, std_err = linregress(x_log, y_log)
                line_x = np.array([x_log.min(), x_log.max()])
                line_y = slope * line_x + intercept
                ax.plot(line_x, line_y, 'r--', linewidth=2, label=f'R² = {r_value**2:.3f}')
                ax.legend(fontsize=11)
                plt.tight_layout()
                self._save_figure('17_scatter_spots_vs_bases')
                plt.close()
    
    def analyze_study_level(self):
        """Analyze samples per study."""
        print("\n=== Analyzing Study-Level Data ===")
        
        study_col = None
        for col in ['study_accession', 'study', 'bioproject']:
            if col in self.df.columns:
                study_col = col
                break
        
        if study_col:
            samples_per_study = self.df[study_col].value_counts()
            
            # Study size distribution
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.hist(samples_per_study.values, bins=50, color='teal', edgecolor='black',
                    alpha=0.7, linewidth=1.2)
            ax.set_xlabel('Number of Samples per Study', fontsize=12, fontweight='bold')
            ax.set_ylabel('Frequency (Number of Studies)', fontsize=12, fontweight='bold')
            ax.set_title('Distribution of Study Sizes', fontsize=14, fontweight='bold', pad=20)
            ax.grid(alpha=0.3, linestyle='--')
            ax.set_yscale('log')
            plt.tight_layout()
            self._save_figure('18_study_size_distribution')
            plt.close()
            
            # Top studies
            top_studies = samples_per_study.head(15)
            
            fig, ax = plt.subplots(figsize=(14, 8))
            colors = sns.color_palette("magma", n_colors=len(top_studies))
            bars = ax.barh(range(len(top_studies)), top_studies.values, color=colors,
                           edgecolor='black', linewidth=1.2)
            ax.set_yticks(range(len(top_studies)))
            ax.set_yticklabels([str(study)[:30] for study in top_studies.index], fontsize=9)
            ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Top 15 Largest Studies by Sample Count', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            for i, val in enumerate(top_studies.values):
                ax.text(val + max(top_studies.values)*0.01, i, f'{val:,}', 
                        va='center', fontsize=9, fontweight='bold')
            ax.invert_yaxis()
            plt.tight_layout()
            self._save_figure('19_top_studies')
            plt.close()
    
    def analyze_library_layout(self):
        """Compare PAIRED vs SINGLE end sequencing."""
        print("\n=== Analyzing Library Layout Comparison ===")
        
        if 'library_layout' in self.df.columns and 'total_spots_numeric' in self.df.columns:
            layout_depth = self.df[['library_layout', 'total_spots_numeric']].dropna()
            layout_depth = layout_depth[layout_depth['total_spots_numeric'] > 0]
            
            if len(layout_depth) > 10:
                layouts = layout_depth['library_layout'].unique()
                data_by_layout = [np.log10(layout_depth[layout_depth['library_layout'] == layout]['total_spots_numeric']) 
                                 for layout in layouts]
                
                # Violin plot
                fig, ax = plt.subplots(figsize=(12, 7))
                violin_parts = ax.violinplot(data_by_layout, positions=range(len(layouts)), 
                                             widths=0.7, showmeans=True, showmedians=True)
                
                colors = sns.color_palette("Set2", n_colors=len(layouts))
                for i, pc in enumerate(violin_parts['bodies']):
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.7)
                
                ax.set_xticks(range(len(layouts)))
                ax.set_xticklabels(layouts, fontsize=11)
                ax.set_ylabel('Log10(Total Spots)', fontsize=12, fontweight='bold')
                ax.set_xlabel('Library Layout', fontsize=12, fontweight='bold')
                ax.set_title('Sequencing Depth Distribution by Library Layout', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('20_library_layout_violin')
                plt.close()
        
        # Boxplot for read length by layout
        if 'library_layout' in self.df.columns and 'avg_read_length_numeric' in self.df.columns:
            layout_length = self.df[['library_layout', 'avg_read_length_numeric']].dropna()
            layout_length = layout_length[layout_length['avg_read_length_numeric'] > 0]
            
            if len(layout_length) > 10:
                layouts = layout_length['library_layout'].unique()
                data_by_layout = [layout_length[layout_length['library_layout'] == layout]['avg_read_length_numeric'] 
                                 for layout in layouts]
                
                fig, ax = plt.subplots(figsize=(12, 7))
                bp = ax.boxplot(data_by_layout, labels=layouts, patch_artist=True, widths=0.6,
                                boxprops=dict(linewidth=1.5, edgecolor='black'),
                                medianprops=dict(color='red', linewidth=2),
                                whiskerprops=dict(linewidth=1.5),
                                capprops=dict(linewidth=1.5))
                
                colors = sns.color_palette("pastel", n_colors=len(layouts))
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                
                ax.set_ylabel('Average Read Length (bp)', fontsize=12, fontweight='bold')
                ax.set_xlabel('Library Layout', fontsize=12, fontweight='bold')
                ax.set_title('Read Length Distribution by Library Layout', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('21_library_layout_boxplot')
                plt.close()
        
        # Boxplot for read length
        if 'library_layout' in self.df.columns and 'avg_read_length_numeric' in self.df.columns:
            layout_length = self.df[['library_layout', 'avg_read_length_numeric']].dropna()
            layout_length = layout_length[layout_length['avg_read_length_numeric'] > 0]
            
            if len(layout_length) > 10:
                layouts = layout_length['library_layout'].unique()
                data_by_layout = [layout_length[layout_length['library_layout'] == layout]['avg_read_length_numeric'] 
                                 for layout in layouts]
                
                fig, ax = plt.subplots(figsize=(12, 7))
                bp = ax.boxplot(data_by_layout, labels=layouts, patch_artist=True, widths=0.6,
                                boxprops=dict(linewidth=1.5, edgecolor='black'),
                                medianprops=dict(color='red', linewidth=2),
                                whiskerprops=dict(linewidth=1.5),
                                capprops=dict(linewidth=1.5))
                
                colors = sns.color_palette("pastel", n_colors=len(layouts))
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                
                ax.set_ylabel('Average Read Length (bp)', fontsize=12, fontweight='bold')
                ax.set_xlabel('Library Layout', fontsize=12, fontweight='bold')
                ax.set_title('Read Length Distribution by Library Layout', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('21_library_layout_boxplot')
                plt.close()
    
    def analyze_patient_demographics_detailed(self):
        """Detailed patient demographics analysis."""
        print("\n=== Analyzing Detailed Patient Demographics ===")
        
        patient_demo_keywords = ['age', 'sex', 'gender', 'race', 'ethnicity', 'population', 
                                 'geographic', 'country', 'nationality']
        patient_demo_cols = []
        for col in self.df.columns:
            col_lower = col.lower()
            if any(keyword in col_lower for keyword in patient_demo_keywords):
                if self.df[col].notna().sum() > 10:
                    patient_demo_cols.append(col)
        
        patient_demo_summary = {}
        for col in patient_demo_cols:
            non_null = self.df[col].notna().sum()
            unique_vals = self.df[col].nunique()
            completeness = (non_null / len(self.df)) * 100
            patient_demo_summary[col] = {
                'non_null': non_null,
                'unique_values': unique_vals,
                'completeness': completeness
            }
        
        if len(patient_demo_summary) > 0:
            demo_df = pd.DataFrame(patient_demo_summary).T
            demo_df = demo_df.sort_values('completeness', ascending=False)
            
            fig, ax = plt.subplots(figsize=(14, 10))
            n_show = min(30, len(demo_df))
            top_demo = demo_df.head(n_show)
            
            colors = plt.cm.RdYlGn(top_demo['completeness'].values / 100)
            bars = ax.barh(range(len(top_demo)), top_demo['completeness'].values,
                           color=colors, edgecolor='black', linewidth=0.8)
            ax.set_yticks(range(len(top_demo)))
            ax.set_yticklabels([col.replace('sample_attributes_', '').replace('_', ' ').title() 
                                for col in top_demo.index], fontsize=9)
            ax.set_xlabel('Completeness (%)', fontsize=12, fontweight='bold')
            ax.set_title('Patient Demographics: Data Completeness', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            ax.set_xlim(0, 100)
            
            for i, (idx, row) in enumerate(top_demo.iterrows()):
                ax.text(row['completeness'] + 1, i, 
                        f"{row['completeness']:.1f}% (n={int(row['non_null']):,})", 
                        va='center', fontsize=7)
            ax.invert_yaxis()
            plt.tight_layout()
            self._save_figure('22_patient_demographics_completeness')
            plt.close()
    
    def analyze_completeness_extended(self):
        """Extended metadata completeness analysis."""
        print("\n=== Extended Completeness Analysis ===")
        
        completeness = (1 - self.df.isnull().sum() / len(self.df)) * 100
        
        # Completeness distribution histogram
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.hist(completeness.values, bins=50, color='teal', edgecolor='black',
                alpha=0.7, linewidth=1.2)
        ax.set_xlabel('Completeness (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Columns', fontsize=12, fontweight='bold')
        ax.set_title('Distribution of Column Completeness', fontsize=14, fontweight='bold', pad=20)
        ax.grid(alpha=0.3, linestyle='--')
        plt.tight_layout()
        self._save_figure('23_completeness_distribution')
        plt.close()
        
        # Completeness categories
        completeness_categories = pd.cut(completeness, bins=[0, 25, 50, 75, 100], 
                                          labels=['0-25%', '25-50%', '50-75%', '75-100%'])
        category_counts = completeness_categories.value_counts().sort_index()
        
        fig, ax = plt.subplots(figsize=(10, 6))
        colors_cat = sns.color_palette("RdYlGn", n_colors=len(category_counts))
        bars = ax.bar(range(len(category_counts)), category_counts.values, color=colors_cat,
                      edgecolor='black', linewidth=1.5, alpha=0.8)
        ax.set_xticks(range(len(category_counts)))
        ax.set_xticklabels(category_counts.index, fontsize=11)
        ax.set_ylabel('Number of Columns', fontsize=12, fontweight='bold')
        ax.set_xlabel('Completeness Range', fontsize=12, fontweight='bold')
        ax.set_title('Metadata Completeness by Category', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        for i, val in enumerate(category_counts.values):
            ax.text(i, val + 2, f'{val}', ha='center', fontsize=11, fontweight='bold')
        plt.tight_layout()
        self._save_figure('24_completeness_categories')
        plt.close()
    
    def analyze_column_categories(self):
        """Analyze column distribution by prefix categories."""
        print("\n=== Analyzing Column Categories ===")
        
        column_categories = {
            'run_': [],
            'experiment_': [],
            'sample_': [],
            'study_': [],
            'submission_': [],
            'organization_': [],
            'other': []
        }
        
        for col in self.df.columns:
            categorized = False
            for prefix in ['run_', 'experiment_', 'sample_', 'study_', 'submission_', 'organization_']:
                if col.startswith(prefix):
                    column_categories[prefix].append(col)
                    categorized = True
                    break
            if not categorized:
                column_categories['other'].append(col)
        
        category_counts = {k: len(v) for k, v in column_categories.items()}
        
        # Bar chart
        fig, ax = plt.subplots(figsize=(12, 7))
        categories = list(category_counts.keys())
        counts = list(category_counts.values())
        colors = sns.color_palette("viridis", n_colors=len(categories))
        bars = ax.bar(range(len(categories)), counts, color=colors, edgecolor='black', 
                      linewidth=1.5, alpha=0.8)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(categories, fontsize=11, rotation=45, ha='right')
        ax.set_ylabel('Number of Columns', fontsize=12, fontweight='bold')
        ax.set_title('Column Distribution by Category Prefix', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        for i, count in enumerate(counts):
            ax.text(i, count + 2, f'{count}', ha='center', fontsize=11, fontweight='bold')
        plt.tight_layout()
        self._save_figure('25_column_categories_bar')
        plt.close()
        
        # Pie chart
        fig, ax = plt.subplots(figsize=(10, 10))
        sizes = [v for v in category_counts.values() if v > 0]
        labels_pie = [f'{k}\n({v} cols)' for k, v in category_counts.items() if v > 0]
        colors_pie = sns.color_palette("Spectral", n_colors=len(sizes))
        wedges, texts, autotexts = ax.pie(sizes, labels=labels_pie, autopct='%1.1f%%', startangle=90,
                                            colors=colors_pie, textprops={'fontsize': 10},
                                            pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(11)
        ax.set_title('Column Category Distribution', fontsize=14, fontweight='bold', pad=20)
        plt.tight_layout()
        self._save_figure('26_column_categories_pie')
        plt.close()
    
    def analyze_technical_parameters(self):
        """Analyze technical sequencing parameters."""
        print("\n=== Analyzing Technical Parameters ===")
        
        # Library layout
        if 'library_layout' in self.df.columns:
            layout_data = self.df['library_layout'].value_counts()
            if len(layout_data) > 0:
                fig, ax = plt.subplots(figsize=(10, 6))
                colors = sns.color_palette("Set3", n_colors=len(layout_data))
                bars = ax.bar(range(len(layout_data)), layout_data.values, color=colors,
                              edgecolor='black', linewidth=1.5, alpha=0.8)
                ax.set_xticks(range(len(layout_data)))
                ax.set_xticklabels(layout_data.index, fontsize=11)
                ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Library Layout Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                for i, val in enumerate(layout_data.values):
                    ax.text(i, val + 50, f'{val:,}', ha='center', fontsize=11, fontweight='bold')
                plt.tight_layout()
                self._save_figure('27_library_layout')
                plt.close()
        
        # Library selection
        if 'library_selection' in self.df.columns:
            selection_data = self.df['library_selection'].value_counts()
            if len(selection_data) > 0:
                fig, ax = plt.subplots(figsize=(12, 7))
                colors = sns.color_palette("viridis", n_colors=len(selection_data))
                bars = ax.barh(range(len(selection_data)), selection_data.values, color=colors,
                               edgecolor='black', linewidth=1.2)
                ax.set_yticks(range(len(selection_data)))
                ax.set_yticklabels(selection_data.index, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Library Selection Method Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                for i, val in enumerate(selection_data.values):
                    ax.text(val + 20, i, f'{val:,}', va='center', fontsize=9, fontweight='bold')
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('28_library_selection')
                plt.close()
        
        # Library source
        if 'library_source' in self.df.columns:
            source_data = self.df['library_source'].value_counts()
            if len(source_data) > 0:
                fig, ax = plt.subplots(figsize=(10, 6))
                colors = sns.color_palette("RdYlGn", n_colors=len(source_data))
                bars = ax.bar(range(len(source_data)), source_data.values, color=colors,
                              edgecolor='black', linewidth=1.5, alpha=0.8)
                ax.set_xticks(range(len(source_data)))
                ax.set_xticklabels(source_data.index, fontsize=11, rotation=45, ha='right')
                ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Library Source Distribution', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                for i, val in enumerate(source_data.values):
                    ax.text(i, val + 50, f'{val:,}', ha='center', fontsize=11, fontweight='bold')
                plt.tight_layout()
                self._save_figure('29_library_source')
                plt.close()
        
        # Platform type pie
        if 'platform_type' in self.df.columns:
            platform_data = self.df['platform_type'].value_counts()
            if len(platform_data) > 0:
                fig, ax = plt.subplots(figsize=(10, 10))
                sizes = platform_data.values
                labels = [f'{plat}\n({cnt:,})' for plat, cnt in zip(platform_data.index, sizes)]
                colors_pie = sns.color_palette("Paired", n_colors=len(platform_data))
                wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90,
                                                    colors=colors_pie, textprops={'fontsize': 10},
                                                    pctdistance=0.85, wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(11)
                ax.set_title('Platform Type Distribution', fontsize=14, fontweight='bold', pad=20)
                plt.tight_layout()
                self._save_figure('30_platform_type_pie')
                plt.close()
    
    def analyze_quality_metrics(self):
        """Analyze quality and QC metrics."""
        print("\n=== Analyzing Quality Metrics ===")
        
        # MBases distribution
        if 'mbases' in self.df.columns:
            mbases_data = pd.to_numeric(self.df['mbases'], errors='coerce').dropna()
            if len(mbases_data) > 0:
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.hist(np.log10(mbases_data + 1), bins=50, color='darkgreen', 
                        edgecolor='black', alpha=0.7, linewidth=1.2)
                ax.set_xlabel('Log10(MBases + 1)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
                ax.set_title('Distribution of Sequencing Yield (MBases)', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('31_mbases_distribution')
                plt.close()
        
        # Core metadata completeness
        core_metadata_cols = ['instrument_model', 'library_strategy', 'library_source', 
                              'library_layout', 'platform_type']
        core_available = [col for col in core_metadata_cols if col in self.df.columns]
        core_completeness = {}
        for col in core_available:
            comp = (1 - self.df[col].isnull().sum() / len(self.df)) * 100
            core_completeness[col] = comp
        
        if len(core_completeness) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            cols = list(core_completeness.keys())
            vals = list(core_completeness.values())
            colors = plt.cm.RdYlGn([v/100 for v in vals])
            bars = ax.bar(range(len(cols)), vals, color=colors, edgecolor='black',
                          linewidth=1.5, alpha=0.8)
            ax.set_xticks(range(len(cols)))
            ax.set_xticklabels([c.replace('_', ' ').title() for c in cols], 
                                fontsize=10, rotation=45, ha='right')
            ax.set_ylabel('Completeness (%)', fontsize=12, fontweight='bold')
            ax.set_title('Core Metadata Completeness', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            ax.set_ylim(0, 100)
            for i, val in enumerate(vals):
                ax.text(i, val + 2, f'{val:.1f}%', ha='center', fontsize=10, fontweight='bold')
            plt.tight_layout()
            self._save_figure('32_core_metadata_completeness')
            plt.close()
    
    def analyze_statistical_summary(self):
        """Generate statistical summary heatmap."""
        print("\n=== Generating Statistical Summary ===")
        
        stat_cols = ['total_spots_numeric', 'total_bases_numeric', 'avg_read_length_numeric']
        available_stats = [col for col in stat_cols if col in self.df.columns and self.df[col].notna().sum() > 0]
        
        if len(available_stats) > 0:
            stats_summary = []
            for col in available_stats:
                data = self.df[col].dropna()
                stats_summary.append({
                    'Metric': col.replace('_numeric', '').replace('_', ' ').title(),
                    'Count': len(data),
                    'Mean': data.mean(),
                    'Median': data.median(),
                    'Std': data.std(),
                    'Min': data.min(),
                    'Max': data.max(),
                    'Q25': data.quantile(0.25),
                    'Q75': data.quantile(0.75)
                })
            
            stats_df = pd.DataFrame(stats_summary)
            
            # Heatmap
            stats_normalized = stats_df.set_index('Metric').iloc[:, 1:]
            stats_normalized_scaled = (stats_normalized - stats_normalized.min()) / (stats_normalized.max() - stats_normalized.min())
            
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.heatmap(stats_normalized_scaled.T, annot=stats_normalized.T, fmt='.2e', 
                        cmap='YlOrRd', cbar_kws={'label': 'Normalized Value'}, ax=ax,
                        linewidths=1, linecolor='white')
            ax.set_title('Statistical Summary Heatmap', fontsize=14, fontweight='bold', pad=20)
            ax.set_xlabel('Metric', fontsize=12, fontweight='bold')
            ax.set_ylabel('Statistic', fontsize=12, fontweight='bold')
            plt.tight_layout()
            self._save_figure('33_statistical_summary_heatmap')
            plt.close()
    
    def create_final_summary_chart(self):
        """Create final dataset summary visualization."""
        print("\n=== Creating Final Summary Chart ===")
        
        summary_stats = {
            'total_samples': len(self.df),
            'total_columns': len(self.df.columns),
            'sample_attributes': len([c for c in self.df.columns if c.startswith('sample_')]),
            'run_attributes': len([c for c in self.df.columns if c.startswith('run_')]),
            'experiment_attributes': len([c for c in self.df.columns if c.startswith('experiment_')]),
            'study_attributes': len([c for c in self.df.columns if c.startswith('study_')])
        }
        
        fig, ax = plt.subplots(figsize=(10, 6))
        categories = list(summary_stats.keys())
        values = list(summary_stats.values())
        colors = sns.color_palette("Set1", n_colors=len(categories))
        bars = ax.bar(range(len(categories)), values, color=colors, edgecolor='black',
                      linewidth=1.5, alpha=0.8)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels([c.replace('_', ' ').title() for c in categories], 
                            fontsize=11, rotation=45, ha='right')
        ax.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax.set_title('Dataset Summary Statistics', fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        for i, val in enumerate(values):
            ax.text(i, val + max(values)*0.02, f'{val:,}', ha='center', fontsize=11, fontweight='bold')
        plt.tight_layout()
        self._save_figure('34_dataset_summary')
        plt.close()
    
    def analyze_patient_attributes(self):
        """Analyze patient-specific attributes comprehensively."""
        print("\n=== Analyzing Patient Attributes ===")
        
        # Get all patient-related columns
        patient_cols = [col for col in self.df.columns if any(x in col.lower() 
                       for x in ['patient', 'donor', 'subject', 'individual'])]
        
        if len(patient_cols) > 0:
            # Calculate completeness for patient columns
            patient_completeness = {}
            for col in patient_cols:
                non_null = self.df[col].notna().sum()
                patient_completeness[col] = non_null
            
            # Sort and take top 40
            patient_completeness = dict(sorted(patient_completeness.items(), 
                                              key=lambda x: x[1], reverse=True)[:40])
            
            if len(patient_completeness) > 0:
                # Bar chart for top 40
                fig, ax = plt.subplots(figsize=(14, 10))
                y_pos = np.arange(len(patient_completeness))
                values = list(patient_completeness.values())
                labels = [col.replace('sample_attributes_', '').replace('_', ' ')[:50] 
                         for col in patient_completeness.keys()]
                
                colors = sns.color_palette("viridis", n_colors=len(values))
                bars = ax.barh(y_pos, values, color=colors, edgecolor='black', linewidth=0.8, alpha=0.85)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(labels, fontsize=9)
                ax.set_xlabel('Number of Non-Null Values', fontsize=12, fontweight='bold')
                ax.set_title('Top 40 Patient Attributes by Data Availability', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, v in enumerate(values):
                    ax.text(v + max(values)*0.01, i, f'{v:,}', va='center', fontsize=8, fontweight='bold')
                
                plt.tight_layout()
                self._save_figure('35_patient_attributes_top40')
                plt.close()
                
                # Histogram of completeness distribution
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.hist(values, bins=30, color='teal', edgecolor='black', alpha=0.7, linewidth=1.2)
                ax.set_xlabel('Number of Non-Null Values', fontsize=12, fontweight='bold')
                ax.set_ylabel('Frequency (Number of Attributes)', fontsize=12, fontweight='bold')
                ax.set_title('Distribution of Patient Attribute Completeness', fontsize=14, fontweight='bold', pad=20)
                ax.grid(alpha=0.3, linestyle='--')
                plt.tight_layout()
                self._save_figure('36_patient_attributes_histogram')
                plt.close()
    
    def analyze_clinical_attributes(self):
        """Analyze clinical and medical attributes."""
        print("\\n=== Analyzing Clinical Attributes ===")
        
        # Find clinical columns
        clinical_keywords = ['disease', 'diagnosis', 'condition', 'treatment', 'therapy', 
                            'medication', 'symptom', 'clinical', 'phenotype', 'stage', 'grade']
        clinical_cols = [col for col in self.df.columns if any(kw in col.lower() for kw in clinical_keywords)]
        
        if len(clinical_cols) > 0:
            # Calculate availability
            clinical_availability = {}
            for col in clinical_cols:
                non_null_pct = (self.df[col].notna().sum() / len(self.df)) * 100
                clinical_availability[col] = non_null_pct
            
            clinical_availability = dict(sorted(clinical_availability.items(), 
                                               key=lambda x: x[1], reverse=True)[:25])
            
            if len(clinical_availability) > 0:
                fig, ax = plt.subplots(figsize=(14, 8))
                y_pos = np.arange(len(clinical_availability))
                values = list(clinical_availability.values())
                labels = [col.replace('sample_attributes_', '').replace('_', ' ')[:50] 
                         for col in clinical_availability.keys()]
                
                colors = sns.color_palette("RdYlGn", n_colors=len(values))
                bars = ax.barh(y_pos, values, color=colors, edgecolor='black', linewidth=0.8, alpha=0.85)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(labels, fontsize=9)
                ax.set_xlabel('Data Availability (%)', fontsize=12, fontweight='bold')
                ax.set_title('Clinical Attributes Data Availability', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, v in enumerate(values):
                    ax.text(v + 1, i, f'{v:.1f}%', va='center', fontsize=8, fontweight='bold')
                
                plt.tight_layout()
                self._save_figure('37_clinical_attributes')
                plt.close()
    
    def analyze_data_collection(self):
        """Analyze data collection and submission metadata."""
        print("\n=== Analyzing Data Collection Metadata ===")
        
        # Dual panel: Centers and Submitters
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Panel 1: Center names
        if 'center_name' in self.df.columns:
            center_counts = self.df['center_name'].value_counts().head(15)
            if len(center_counts) > 0:
                colors1 = sns.color_palette("tab20", n_colors=len(center_counts))
                bars1 = ax1.barh(range(len(center_counts)), center_counts.values, 
                                color=colors1, edgecolor='black', linewidth=1.2, alpha=0.8)
                ax1.set_yticks(range(len(center_counts)))
                ax1.set_yticklabels([str(c)[:40] for c in center_counts.index], fontsize=9)
                ax1.set_xlabel('Number of Samples', fontsize=11, fontweight='bold')
                ax1.set_title('Top 15 Data Collection Centers', fontsize=12, fontweight='bold')
                ax1.grid(axis='x', alpha=0.3, linestyle='--')
                ax1.invert_yaxis()
                
                for i, v in enumerate(center_counts.values):
                    ax1.text(v + max(center_counts.values)*0.01, i, f'{v:,}', 
                            va='center', fontsize=8, fontweight='bold')
        
        # Panel 2: Submission dates
        if 'run_submission_date' in self.df.columns:
            self.df['submission_year'] = pd.to_datetime(self.df['run_submission_date'], errors='coerce').dt.year
            year_counts = self.df['submission_year'].value_counts().sort_index()
            
            if len(year_counts) > 0:
                colors2 = sns.color_palette("coolwarm", n_colors=len(year_counts))
                bars2 = ax2.bar(year_counts.index, year_counts.values, 
                               color=colors2, edgecolor='black', linewidth=1.2, alpha=0.8)
                ax2.set_xlabel('Submission Year', fontsize=11, fontweight='bold')
                ax2.set_ylabel('Number of Samples', fontsize=11, fontweight='bold')
                ax2.set_title('Data Submission Timeline', fontsize=12, fontweight='bold')
                ax2.grid(axis='y', alpha=0.3, linestyle='--')
                plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        self._save_figure('38_data_collection_overview')
        plt.close()
    
    def analyze_sample_collection(self):
        """Analyze sample collection metadata."""
        print("\n=== Analyzing Sample Collection ===")
        
        # Look for collection-related columns
        collection_cols = [col for col in self.df.columns if any(x in col.lower() 
                          for x in ['collection', 'sampling', 'specimen', 'biopsy'])]
        
        if len(collection_cols) > 0:
            # Collection site/tissue analysis
            tissue_col = None
            for col in ['sample_attributes_tissue', 'tissue', 'sample_attributes_source_name']:
                if col in self.df.columns:
                    tissue_col = col
                    break
            
            if tissue_col:
                tissue_counts = self.df[tissue_col].value_counts().head(20)
                
                fig, ax = plt.subplots(figsize=(14, 8))
                colors = sns.color_palette("Spectral", n_colors=len(tissue_counts))
                bars = ax.barh(range(len(tissue_counts)), tissue_counts.values,
                              color=colors, edgecolor='black', linewidth=1.2, alpha=0.85)
                ax.set_yticks(range(len(tissue_counts)))
                ax.set_yticklabels([str(t)[:50] for t in tissue_counts.index], fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sample Collection: Top 20 Tissue/Source Types', fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                ax.invert_yaxis()
                
                for i, v in enumerate(tissue_counts.values):
                    ax.text(v + max(tissue_counts.values)*0.01, i, f'{v:,}', 
                           va='center', fontsize=9, fontweight='bold')
                
                plt.tight_layout()
                self._save_figure('39_sample_collection_tissues')
                plt.close()
    
    def analyze_missing_data_patterns(self):
        """Analyze missing data patterns across columns."""
        print("\n=== Analyzing Missing Data Patterns ===")
        
        # Calculate missingness
        missing_pct = (self.df.isnull().sum() / len(self.df) * 100).sort_values(ascending=False)
        missing_top = missing_pct[missing_pct > 50].head(30)
        
        if len(missing_top) > 0:
            fig, ax = plt.subplots(figsize=(14, 10))
            y_pos = np.arange(len(missing_top))
            values = missing_top.values
            labels = [col.replace('sample_attributes_', '').replace('_', ' ')[:50] 
                     for col in missing_top.index]
            
            colors = sns.color_palette("YlOrRd", n_colors=len(values))
            bars = ax.barh(y_pos, values, color=colors, edgecolor='black', linewidth=0.8, alpha=0.85)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(labels, fontsize=9)
            ax.set_xlabel('Missing Data (%)', fontsize=12, fontweight='bold')
            ax.set_title('Top 30 Columns with Highest Missing Data', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            ax.invert_yaxis()
            
            for i, v in enumerate(values):
                ax.text(v + 1, i, f'{v:.1f}%', va='center', fontsize=8, fontweight='bold')
            
            plt.tight_layout()
            self._save_figure('40_missing_data_heatmap')
            plt.close()
    
    def analyze_complete_patient_profile(self):
        """Generate comprehensive patient profile summary."""
        print("\n=== Creating Complete Patient Profile ===")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
        
        # Panel 1: Patient attribute categories
        patient_categories = {
            'Demographics': 0,
            'Clinical': 0,
            'Sample Info': 0,
            'Treatment': 0,
            'Outcome': 0
        }
        
        for col in self.df.columns:
            col_lower = col.lower()
            if any(x in col_lower for x in ['age', 'sex', 'gender', 'race', 'ethnicity']):
                patient_categories['Demographics'] += 1
            elif any(x in col_lower for x in ['disease', 'diagnosis', 'condition', 'phenotype', 'stage']):
                patient_categories['Clinical'] += 1
            elif any(x in col_lower for x in ['tissue', 'specimen', 'biopsy', 'collection']):
                patient_categories['Sample Info'] += 1
            elif any(x in col_lower for x in ['treatment', 'therapy', 'drug', 'medication']):
                patient_categories['Treatment'] += 1
            elif any(x in col_lower for x in ['outcome', 'survival', 'response', 'prognosis']):
                patient_categories['Outcome'] += 1
        
        # Remove empty categories
        patient_categories = {k: v for k, v in patient_categories.items() if v > 0}
        
        if len(patient_categories) > 0:
            colors1 = sns.color_palette("Set3", n_colors=len(patient_categories))
            bars1 = ax1.bar(range(len(patient_categories)), patient_categories.values(),
                           color=colors1, edgecolor='black', linewidth=1.5, alpha=0.8)
            ax1.set_xticks(range(len(patient_categories)))
            ax1.set_xticklabels(patient_categories.keys(), fontsize=11, fontweight='bold')
            ax1.set_ylabel('Number of Attributes', fontsize=12, fontweight='bold')
            ax1.set_title('Patient Data Categories', fontsize=13, fontweight='bold')
            ax1.grid(axis='y', alpha=0.3, linestyle='--')
            
            for i, (cat, val) in enumerate(patient_categories.items()):
                ax1.text(i, val + max(patient_categories.values())*0.02, str(val),
                        ha='center', fontsize=11, fontweight='bold')
        
        # Panel 2: Data completeness summary
        completeness_summary = {
            'Highly Complete (>80%)': 0,
            'Moderately Complete (50-80%)': 0,
            'Sparsely Complete (20-50%)': 0,
            'Mostly Missing (<20%)': 0
        }
        
        for col in self.df.columns:
            completeness = (self.df[col].notna().sum() / len(self.df)) * 100
            if completeness > 80:
                completeness_summary['Highly Complete (>80%)'] += 1
            elif completeness > 50:
                completeness_summary['Moderately Complete (50-80%)'] += 1
            elif completeness > 20:
                completeness_summary['Sparsely Complete (20-50%)'] += 1
            else:
                completeness_summary['Mostly Missing (<20%)'] += 1
        
        colors2 = ['#2ecc71', '#f39c12', '#e67e22', '#e74c3c']
        wedges, texts, autotexts = ax2.pie(completeness_summary.values(), 
                                            labels=[f'{k}\\n({v})' for k, v in completeness_summary.items()],
                                            autopct='%1.1f%%', startangle=90, colors=colors2,
                                            textprops={'fontsize': 10}, pctdistance=0.85,
                                            wedgeprops=dict(width=0.6, edgecolor='white', linewidth=2))
        
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(11)
        
        ax2.set_title('Overall Data Completeness', fontsize=13, fontweight='bold')
        
        plt.tight_layout()
        self._save_figure('41_complete_patient_profile')
        plt.close()
    
    def analyze_categorical_patients(self):
        """Multi-panel analysis of categorical patient variables."""
        print("\n=== Analyzing Categorical Patient Variables ===")
        
        # Find categorical patient columns
        categorical_cols = []
        for col in self.df.columns:
            if any(x in col.lower() for x in ['patient', 'sample_attributes']) and self.df[col].dtype == 'object':
                unique_count = self.df[col].nunique()
                if 2 <= unique_count <= 20:  # Reasonable number of categories
                    categorical_cols.append(col)
        
        if len(categorical_cols) >= 6:
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            axes = axes.flatten()
            
            for idx, col in enumerate(categorical_cols[:6]):
                value_counts = self.df[col].value_counts().head(10)
                
                colors = sns.color_palette("husl", n_colors=len(value_counts))
                bars = axes[idx].bar(range(len(value_counts)), value_counts.values,
                                    color=colors, edgecolor='black', linewidth=1.2, alpha=0.8)
                axes[idx].set_xticks(range(len(value_counts)))
                axes[idx].set_xticklabels([str(v)[:20] for v in value_counts.index], 
                                         rotation=45, ha='right', fontsize=8)
                axes[idx].set_ylabel('Count', fontsize=10, fontweight='bold')
                axes[idx].set_title(col.replace('sample_attributes_', '').replace('_', ' ')[:40],
                                   fontsize=11, fontweight='bold')
                axes[idx].grid(axis='y', alpha=0.3, linestyle='--')
                
                for i, v in enumerate(value_counts.values):
                    axes[idx].text(i, v + max(value_counts.values)*0.02, f'{v:,}',
                                  ha='center', fontsize=8, fontweight='bold')
            
            plt.suptitle('Categorical Patient Variables Distribution', fontsize=16, fontweight='bold', y=0.995)
            plt.tight_layout()
            self._save_figure('42_categorical_patients_multipanel')
            plt.close()
    
    def analyze_sample_level_summary(self):
        """Analyze data at sample level (deduplicated patients)."""
        print("\n=== Sample-Level Analysis (Unique Patients) ===")
        
        # Create sample-level aggregated dataframe
        # For each sample, take the first run and aggregate numeric data
        sample_cols = ['sample_accession', 'sample_sex', 'sample_age', 'sample_tissue', 
                      'sample_disease', 'organism', 'study_accession']
        
        # Count runs per sample
        runs_per_sample = self.df.groupby('sample_accession').size().reset_index(name='run_count')
        
        # Get unique samples
        df_samples = self.df.drop_duplicates(subset='sample_accession').copy()
        df_samples = df_samples.merge(runs_per_sample, on='sample_accession', how='left')
        
        print(f"Total runs: {self.n_runs}")
        print(f"Unique samples: {self.n_samples}")
        print(f"Samples with multiple runs: {(df_samples['run_count'] > 1).sum()}")
        
        # Plot 1: Distribution of runs per sample
        fig, ax = plt.subplots(figsize=(12, 6))
        run_dist = df_samples['run_count'].value_counts().sort_index()
        colors = sns.color_palette("coolwarm", n_colors=len(run_dist))
        bars = ax.bar(run_dist.index, run_dist.values, color=colors, 
                     edgecolor='black', linewidth=1.5, alpha=0.8)
        ax.set_xlabel('Number of Runs per Sample', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
        ax.set_title('Distribution of Sequencing Runs per Sample (Patient)', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        for bar, (idx, val) in zip(bars, run_dist.items()):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:,}\n({val/len(df_samples)*100:.1f}%)',
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        plt.tight_layout()
        self._save_figure('43_runs_per_sample_distribution')
        plt.close()
        
        # Plot 2: Comparison - Runs vs Unique Samples
        fig, ax = plt.subplots(figsize=(10, 6))
        categories = ['Total SRA Runs', 'Unique Samples\n(Patients)', 
                     'Samples with\nMultiple Runs', 'Samples with\nSingle Run']
        values = [self.n_runs, self.n_samples, 
                 (df_samples['run_count'] > 1).sum(),
                 (df_samples['run_count'] == 1).sum()]
        colors = ['#e74c3c', '#3498db', '#f39c12', '#2ecc71']
        
        bars = ax.bar(range(len(categories)), values, color=colors, 
                     edgecolor='black', linewidth=2, alpha=0.85)
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(categories, fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax.set_title('SRA Runs vs Unique Samples (Patient-Level Deduplication)', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        for i, (bar, val) in enumerate(zip(bars, values)):
            height = bar.get_height()
            pct = val / values[0] * 100 if i > 0 else 100
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:,}\n({pct:.1f}%)',
                   ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        self._save_figure('44_runs_vs_samples_comparison')
        plt.close()
        
        # Plot 3: Top samples with most runs
        top_samples = df_samples.nlargest(15, 'run_count')[['sample_accession', 'run_count']]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        colors = sns.color_palette("rocket", n_colors=len(top_samples))
        bars = ax.barh(range(len(top_samples)), top_samples['run_count'].values, 
                      color=colors, edgecolor='black', linewidth=1.5)
        ax.set_yticks(range(len(top_samples)))
        ax.set_yticklabels(top_samples['sample_accession'].values, fontsize=10)
        ax.set_xlabel('Number of Runs', fontsize=12, fontweight='bold')
        ax.set_title('Top 15 Samples with Most Sequencing Runs', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        
        for i, val in enumerate(top_samples['run_count'].values):
            ax.text(val + 0.2, i, f'{val}', va='center', fontsize=10, fontweight='bold')
        
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('45_top_samples_by_runs')
        plt.close()
        
        # Store sample-level dataframe for potential further analysis
        self.df_samples = df_samples
        print(f"✓ Generated 3 sample-level visualizations")
    
    def analyze_comprehensive_phenotypes(self):
        """Comprehensive analysis of all phenotype-related columns."""
        print("\n=== Comprehensive Phenotype Analysis ===")
        
        # Define phenotype-related keywords
        phenotype_keywords = ['phenotype', 'disease', 'condition', 'diagnosis', 'clinical', 
                             'trait', 'sex', 'age', 'tissue', 'body_site', 'health', 
                             'symptom', 'complication', 'affected', 'treatment']
        
        # Find all phenotype columns
        phenotype_cols = [col for col in self.df.columns if any(
            keyword in col.lower() for keyword in phenotype_keywords)]
        
        # Filter to columns with data
        phenotype_data = {}
        for col in phenotype_cols:
            non_null = self.df[col].notna().sum()
            if non_null > 10:  # At least 10 samples
                phenotype_data[col] = {
                    'count': non_null,
                    'pct': non_null / len(self.df) * 100,
                    'unique': self.df[col].nunique()
                }
        
        print(f"Found {len(phenotype_data)} phenotype columns with sufficient data")
        
        # Plot 1: Phenotype data completeness - ALL columns
        fig, ax = plt.subplots(figsize=(14, 12))
        sorted_cols = sorted(phenotype_data.items(), key=lambda x: x[1]['count'], reverse=True)
        
        col_names = [col.replace('sample_', '').replace('exp_attr_', '') for col, _ in sorted_cols]
        counts = [data['count'] for _, data in sorted_cols]
        
        colors = sns.color_palette("viridis", n_colors=len(sorted_cols))
        bars = ax.barh(range(len(col_names)), counts, color=colors, 
                      edgecolor='black', linewidth=1.2, alpha=0.85)
        
        ax.set_yticks(range(len(col_names)))
        ax.set_yticklabels(col_names, fontsize=9)
        ax.set_xlabel('Number of Samples with Data', fontsize=12, fontweight='bold')
        ax.set_title('Phenotype Data Availability - All Columns', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        
        for i, (bar, count) in enumerate(zip(bars, counts)):
            pct = count / len(self.df) * 100
            ax.text(count + 20, i, f'{count:,} ({pct:.1f}%)', 
                   va='center', fontsize=8, fontweight='bold')
        
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('46_phenotype_data_completeness')
        plt.close()
        
        # Plot 2: Sex distribution - ALL values
        if 'sample_sex' in self.df.columns:
            sex_data = self.df['sample_sex'].value_counts()
            
            fig, ax = plt.subplots(figsize=(12, 7))
            colors_sex = sns.color_palette("Set2", n_colors=len(sex_data))
            bars = ax.bar(range(len(sex_data)), sex_data.values, color=colors_sex,
                         edgecolor='black', linewidth=1.5, alpha=0.85)
            
            ax.set_xticks(range(len(sex_data)))
            ax.set_xticklabels(sex_data.index, rotation=45, ha='right', fontsize=11)
            ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Sample Sex Distribution - All Categories', fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='y', alpha=0.3, linestyle='--')
            
            for bar, val in zip(bars, sex_data.values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:,}\n({val/len(self.df)*100:.1f}%)',
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            plt.tight_layout()
            self._save_figure('47_sex_distribution_all')
            plt.close()
        
        # Plot 3: Age distribution
        if 'sample_age' in self.df.columns:
            age_data = self.df['sample_age'].dropna()
            
            # Try to extract numeric ages
            def extract_age(val):
                if pd.isna(val):
                    return None
                try:
                    import re
                    match = re.search(r'\d+\.?\d*', str(val))
                    if match:
                        return float(match.group())
                except:
                    pass
                return None
            
            ages_numeric = age_data.apply(extract_age).dropna()
            
            if len(ages_numeric) > 20:
                fig, ax = plt.subplots(figsize=(12, 7))
                ax.hist(ages_numeric, bins=50, color='steelblue', edgecolor='black', 
                       linewidth=1.2, alpha=0.8)
                ax.set_xlabel('Age (years)', fontsize=12, fontweight='bold')
                ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title(f'Age Distribution (n={len(ages_numeric)} samples)', 
                            fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='y', alpha=0.3, linestyle='--')
                
                mean_age = ages_numeric.mean()
                median_age = ages_numeric.median()
                ax.axvline(mean_age, color='red', linestyle='--', linewidth=2, 
                          label=f'Mean: {mean_age:.1f}')
                ax.axvline(median_age, color='green', linestyle='--', linewidth=2,
                          label=f'Median: {median_age:.1f}')
                ax.legend(fontsize=11)
                
                plt.tight_layout()
                self._save_figure('48_age_distribution_histogram')
                plt.close()
        
        # Plot 4: Disease - sample_disease (ALL values)
        if 'sample_disease' in self.df.columns:
            disease_data = self.df['sample_disease'].value_counts()
            
            fig, ax = plt.subplots(figsize=(14, max(10, len(disease_data) * 0.4)))
            colors = sns.color_palette("Reds_r", n_colors=len(disease_data))
            bars = ax.barh(range(len(disease_data)), disease_data.values, color=colors,
                          edgecolor='black', linewidth=1.2)
            
            ax.set_yticks(range(len(disease_data)))
            ax.set_yticklabels(disease_data.index, fontsize=10)
            ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Sample Disease Distribution - All Diseases', 
                        fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            
            for i, val in enumerate(disease_data.values):
                ax.text(val + max(disease_data.values)*0.01, i, f'{val:,}',
                       va='center', fontsize=9, fontweight='bold')
            
            ax.invert_yaxis()
            plt.tight_layout()
            self._save_figure('49_disease_distribution_all')
            plt.close()
        
        # Plot 5: Tissue distribution - ALL values
        if 'sample_tissue' in self.df.columns:
            tissue_data = self.df['sample_tissue'].value_counts()
            
            fig, ax = plt.subplots(figsize=(14, max(10, len(tissue_data) * 0.3)))
            colors_tissue = sns.color_palette("Spectral", n_colors=len(tissue_data))
            bars = ax.barh(range(len(tissue_data)), tissue_data.values, 
                          color=colors_tissue, edgecolor='black', linewidth=1.5)
            
            ax.set_yticks(range(len(tissue_data)))
            ax.set_yticklabels(tissue_data.index, fontsize=10)
            ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
            ax.set_title('Sample Tissue Distribution - All Tissues', 
                        fontsize=14, fontweight='bold', pad=20)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            
            for i, val in enumerate(tissue_data.values):
                ax.text(val + max(tissue_data.values)*0.01, i, f'{val:,}', 
                       va='center', fontsize=9, fontweight='bold')
            
            ax.invert_yaxis()
            plt.tight_layout()
            self._save_figure('50_tissue_distribution_all')
            plt.close()
        
        # Plot 6: Treatment distribution - ALL values
        if 'sample_treatment' in self.df.columns:
            treatment_data = self.df['sample_treatment'].value_counts()
            
            if len(treatment_data) > 0:
                fig, ax = plt.subplots(figsize=(14, max(10, len(treatment_data) * 0.3)))
                colors_treatment = sns.color_palette("coolwarm", n_colors=len(treatment_data))
                bars = ax.barh(range(len(treatment_data)), treatment_data.values,
                              color=colors_treatment, edgecolor='black', linewidth=1.5)
                
                ax.set_yticks(range(len(treatment_data)))
                labels = [str(t)[:80] + '...' if len(str(t)) > 80 else str(t) 
                         for t in treatment_data.index]
                ax.set_yticklabels(labels, fontsize=9)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sample Treatment Distribution - All Treatments', 
                            fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, val in enumerate(treatment_data.values):
                    ax.text(val + max(treatment_data.values)*0.01, i, f'{val}', 
                           va='center', fontsize=9, fontweight='bold')
                
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('51_treatment_distribution_all')
                plt.close()
        
        # Plot 7: Condition distribution - separate graph
        if 'sample_condition' in self.df.columns:
            condition_data = self.df['sample_condition'].value_counts()
            
            if len(condition_data) > 0:
                fig, ax = plt.subplots(figsize=(12, max(8, len(condition_data) * 0.5)))
                colors = sns.color_palette("Blues_r", n_colors=len(condition_data))
                bars = ax.barh(range(len(condition_data)), condition_data.values, 
                              color=colors, edgecolor='black', linewidth=1.2)
                
                ax.set_yticks(range(len(condition_data)))
                ax.set_yticklabels(condition_data.index, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sample Condition Distribution', 
                            fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, val in enumerate(condition_data.values):
                    ax.text(val + max(condition_data.values)*0.02, i, f'{val}',
                           va='center', fontsize=10, fontweight='bold')
                
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('52_condition_distribution_all')
                plt.close()
        
        # Plot 8: Diagnosis distribution - separate graph
        if 'sample_diagnosis' in self.df.columns:
            diagnosis_data = self.df['sample_diagnosis'].value_counts()
            
            if len(diagnosis_data) > 0:
                fig, ax = plt.subplots(figsize=(12, max(8, len(diagnosis_data) * 0.5)))
                colors = sns.color_palette("Greens_r", n_colors=len(diagnosis_data))
                bars = ax.barh(range(len(diagnosis_data)), diagnosis_data.values,
                              color=colors, edgecolor='black', linewidth=1.2)
                
                ax.set_yticks(range(len(diagnosis_data)))
                labels = [str(t)[:60] + '...' if len(str(t)) > 60 else str(t) 
                         for t in diagnosis_data.index]
                ax.set_yticklabels(labels, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sample Diagnosis Distribution', 
                            fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, val in enumerate(diagnosis_data.values):
                    ax.text(val + max(diagnosis_data.values)*0.02, i, f'{val}',
                           va='center', fontsize=10, fontweight='bold')
                
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('53_diagnosis_distribution_all')
                plt.close()
        
        # Plot 9: Disease state - separate graph
        if 'sample_disease_state' in self.df.columns:
            disease_state_data = self.df['sample_disease_state'].value_counts()
            
            if len(disease_state_data) > 0:
                fig, ax = plt.subplots(figsize=(12, max(8, len(disease_state_data) * 0.5)))
                colors = sns.color_palette("Oranges_r", n_colors=len(disease_state_data))
                bars = ax.barh(range(len(disease_state_data)), disease_state_data.values,
                              color=colors, edgecolor='black', linewidth=1.2)
                
                ax.set_yticks(range(len(disease_state_data)))
                labels = [str(t)[:60] + '...' if len(str(t)) > 60 else str(t) 
                         for t in disease_state_data.index]
                ax.set_yticklabels(labels, fontsize=10)
                ax.set_xlabel('Number of Samples', fontsize=12, fontweight='bold')
                ax.set_title('Sample Disease State Distribution', 
                            fontsize=14, fontweight='bold', pad=20)
                ax.grid(axis='x', alpha=0.3, linestyle='--')
                
                for i, val in enumerate(disease_state_data.values):
                    ax.text(val + max(disease_state_data.values)*0.02, i, f'{val}',
                           va='center', fontsize=10, fontweight='bold')
                
                ax.invert_yaxis()
                plt.tight_layout()
                self._save_figure('54_disease_state_distribution_all')
                plt.close()
        
        # Plot 10: Phenotype diversity (unique values per column)
        fig, ax = plt.subplots(figsize=(14, 12))
        sorted_diversity = sorted(phenotype_data.items(), 
                                 key=lambda x: x[1]['unique'], reverse=True)
        
        col_names_div = [col.replace('sample_', '').replace('exp_attr_', '') 
                        for col, _ in sorted_diversity]
        unique_vals = [data['unique'] for _, data in sorted_diversity]
        
        colors_div = sns.color_palette("plasma", n_colors=len(sorted_diversity))
        bars = ax.barh(range(len(col_names_div)), unique_vals, color=colors_div,
                      edgecolor='black', linewidth=1.2)
        
        ax.set_yticks(range(len(col_names_div)))
        ax.set_yticklabels(col_names_div, fontsize=9)
        ax.set_xlabel('Number of Unique Values', fontsize=12, fontweight='bold')
        ax.set_title('Phenotype Diversity - Unique Values per Column', 
                    fontsize=14, fontweight='bold', pad=20)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        
        for i, val in enumerate(unique_vals):
            ax.text(val + max(unique_vals)*0.02, i, f'{val}', 
                   va='center', fontsize=8, fontweight='bold')
        
        ax.invert_yaxis()
        plt.tight_layout()
        self._save_figure('55_phenotype_diversity_all')
        plt.close()
        
        print(f"✓ Generated {10} comprehensive phenotype visualizations")
    
    def run_full_analysis(self):
        """Run complete analysis pipeline with all visualizations."""
        print("\n" + "="*60)
        print("SRA METADATA COMPREHENSIVE ANALYSIS")
        print("Generating 30+ visualizations...")
        print("="*60)
        
        # Core data loading
        self.load_data()
        self.filter_illumina()
        
        # Sequencing analysis (6 plots)
        self.analyze_instruments()  # 2 plots
        self.analyze_library_strategy()  # 2 plots
        self.analyze_sequencing_depth()  # 1 plot
        self.analyze_demographics()  # 1 plot
        
        # Sample analysis (4 plots)
        self.analyze_diseases()  # 2 plots
        self.analyze_tissues()  # 2 plots
        
        # Metadata quality (4 plots)
        self.analyze_metadata_completeness()  # 1 plot
        self.analyze_completeness_extended()  # 2 plots
        self.analyze_column_categories()  # 2 plots
        
        # Organizations and temporal (4 plots)
        self.analyze_organizations()  # 1 plot
        self.analyze_temporal()  # 2 plots
        self.analyze_correlations()  # 1 plot
        
        # Advanced analysis (6 plots)
        self.analyze_scatter_plots()  # 2 plots
        self.analyze_study_level()  # 2 plots
        self.analyze_library_layout()  # 2 plots
        
        # Technical parameters (5 plots)
        self.analyze_technical_parameters()  # 4 plots
        self.analyze_patient_demographics_detailed()  # 1 plot
        
        # Quality and summary (4 plots)
        self.analyze_quality_metrics()  # 2 plots
        self.analyze_statistical_summary()  # 1 plot
        self.create_final_summary_chart()  # 1 plot
        
        # Sample-level analysis (3 plots) - NEW!
        self.analyze_sample_level_summary()  # 3 plots
        
        # Comprehensive phenotype analysis (8 plots) - NEW!
        self.analyze_comprehensive_phenotypes()  # 10 plots
        
        # Export results
        self.generate_summary_report()
        
        # Count generated figures
        num_figures = len(list(self.figures_dir.glob('*.png')))
        
        print("\n" + "="*60)
        print("✅ ANALYSIS COMPLETE!")
        print(f"📊 Generated {num_figures} visualizations")
        print(f"📁 Results: {self.output_dir}")
        print(f"🖼️  Figures: {self.figures_dir}")
        print("="*60 + "\n")


def main():
    """Main entry point for the analysis script."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Comprehensive SRA Metadata Analysis for Illumina Platform',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sra_metadata_analysis.py
  python sra_metadata_analysis.py --input data/metadata/sra_metadata.csv
  python sra_metadata_analysis.py --output results/
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        default='data/metadata/sra_metadata_complete_20251117_085704.csv',
        help='Path to SRA metadata CSV file'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='data/analysis_results',
        help='Output directory for results and figures'
    )
    
    args = parser.parse_args()
    
    # Run analysis
    analyzer = SRAMetadataAnalyzer(args.input, args.output)
    analyzer.run_full_analysis()


if __name__ == '__main__':
    main()
