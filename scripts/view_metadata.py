#!/usr/bin/env python3
"""
Quick viewer for downloaded metadata.
Shows summary statistics and key patient data.
"""

import sys
import pandas as pd
from pathlib import Path

def view_metadata(csv_file: str):
    """Display metadata summary."""
    
    df = pd.read_csv(csv_file)
    
    print("="*70)
    print(f"METADATA SUMMARY: {Path(csv_file).name}")
    print("="*70)
    print(f"Total samples: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    
    # Get patient columns
    sample_cols = [c for c in df.columns if c.startswith('sample_')]
    
    print(f"\n{'─'*70}")
    print("KEY PATIENT DATA:")
    print(f"{'─'*70}")
    
    # Demographics
    if 'sample_sex' in df.columns:
        print("\n📊 SEX/GENDER:")
        counts = df['sample_sex'].value_counts()
        for val, count in counts.items():
            print(f"  {val}: {count} ({count/len(df)*100:.1f}%)")
        missing = df['sample_sex'].isna().sum()
        if missing > 0:
            print(f"  Missing: {missing} ({missing/len(df)*100:.1f}%)")
    
    # Age
    age_cols = [c for c in sample_cols if 'age' in c.lower()]
    if age_cols:
        print(f"\n📅 AGE DATA:")
        for col in age_cols:
            non_null = df[col].dropna()
            if len(non_null) > 0:
                print(f"  {col.replace('sample_', '')}: {len(non_null)} samples have this data")
                print(f"    Examples: {', '.join(map(str, non_null.head(3).values))}")
    
    # Disease
    disease_cols = [c for c in sample_cols if any(x in c.lower() for x in ['disease', 'phenotype', 'condition'])]
    if disease_cols:
        print(f"\n🏥 DISEASE/PHENOTYPE:")
        for col in disease_cols[:3]:
            non_null = df[col].dropna()
            if len(non_null) > 0:
                print(f"\n  {col.replace('sample_', '')}:")
                counts = df[col].value_counts().head(5)
                for val, count in counts.items():
                    val_str = str(val)[:50] + "..." if len(str(val)) > 50 else str(val)
                    print(f"    • {val_str}: {count}")
    
    # Tissue/Body site
    tissue_cols = [c for c in sample_cols if any(x in c.lower() for x in ['tissue', 'body_site', 'cell_type', 'source_name'])]
    if tissue_cols:
        print(f"\n🔬 TISSUE/CELL TYPE:")
        for col in tissue_cols[:3]:
            non_null = df[col].dropna()
            if len(non_null) > 0:
                print(f"\n  {col.replace('sample_', '')}:")
                counts = df[col].value_counts().head(5)
                for val, count in counts.items():
                    val_str = str(val)[:50] + "..." if len(str(val)) > 50 else str(val)
                    print(f"    • {val_str}: {count}")
    
    # Sequencing platform
    if 'platform_type' in df.columns:
        print(f"\n🧬 SEQUENCING PLATFORM:")
        counts = df['platform_type'].value_counts()
        for val, count in counts.items():
            print(f"  {val}: {count}")
    
    if 'library_strategy' in df.columns:
        print(f"\n📚 LIBRARY STRATEGY:")
        counts = df['library_strategy'].value_counts()
        for val, count in counts.items():
            print(f"  {val}: {count}")
    
    # Data completeness
    print(f"\n{'─'*70}")
    print("DATA COMPLETENESS (top patient attributes):")
    print(f"{'─'*70}")
    
    key_cols = [c for c in sample_cols if any(x in c.lower() for x in 
                ['sex', 'age', 'disease', 'tissue', 'body_site', 'cell_type', 
                 'phenotype', 'genotype', 'treatment'])]
    
    if key_cols:
        completeness_data = []
        for col in key_cols:
            pct = (df[col].notna().sum() / len(df)) * 100
            completeness_data.append((col.replace('sample_', ''), pct, df[col].notna().sum()))
        
        # Sort by completeness
        completeness_data.sort(key=lambda x: x[1], reverse=True)
        
        for col_name, pct, count in completeness_data[:15]:
            bar_length = int(pct / 5)
            bar = "█" * bar_length + "░" * (20 - bar_length)
            print(f"  {col_name:30s} {bar} {pct:5.1f}% ({count}/{len(df)})")
    
    print(f"\n{'─'*70}")
    print(f"Available columns: {len(sample_cols)} patient/sample attributes")
    print(f"To see all columns: df.columns.tolist()")
    print(f"To load in Python: df = pd.read_csv('{csv_file}')")
    print("="*70)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Find most recent metadata file
        metadata_dir = Path(__file__).parent.parent / "data" / "metadata"
        csv_files = list(metadata_dir.glob("sra_metadata_complete_*.csv"))
        
        if not csv_files:
            print("No metadata files found!")
            print("Run: python scripts/download_all_metadata.py first")
            sys.exit(1)
        
        # Get most recent
        csv_file = max(csv_files, key=lambda p: p.stat().st_mtime)
        print(f"Using most recent file: {csv_file.name}\n")
    else:
        csv_file = sys.argv[1]
    
    view_metadata(csv_file)
