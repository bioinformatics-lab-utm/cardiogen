#!/usr/bin/env python3
"""
Find samples with most complete metadata.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path


def find_complete_samples(csv_file: str, top_n: int = 5, 
                         min_completeness: float = None,
                         required_fields: list = None):
    """
    Find samples with most complete metadata.
    
    Args:
        csv_file: Path to metadata CSV
        top_n: Number of top samples to return
        min_completeness: Minimum % of fields filled (0-100)
        required_fields: List of required field names (without 'sample_' prefix)
    """
    
    df = pd.read_csv(csv_file, low_memory=False)
    
    # Get sample columns
    sample_cols = [c for c in df.columns if c.startswith('sample_')]
    
    print(f"Total samples: {len(df)}")
    print(f"Total patient attributes: {len(sample_cols)}\n")
    
    # Calculate completeness
    df['completeness'] = df[sample_cols].notna().sum(axis=1)
    df['completeness_pct'] = (df['completeness'] / len(sample_cols)) * 100
    
    # Filter by required fields if specified
    if required_fields:
        print(f"Filtering by required fields: {', '.join(required_fields)}")
        for field in required_fields:
            col_name = f'sample_{field}' if not field.startswith('sample_') else field
            if col_name in df.columns:
                df = df[df[col_name].notna()]
        print(f"Samples after filtering: {len(df)}\n")
    
    # Filter by minimum completeness
    if min_completeness:
        df = df[df['completeness_pct'] >= min_completeness]
        print(f"Samples with >={min_completeness}% completeness: {len(df)}\n")
    
    if len(df) == 0:
        print("No samples match the criteria!")
        return None
    
    # Get top N
    top_samples = df.nlargest(min(top_n, len(df)), 'completeness')
    
    print("="*80)
    print(f"TOP {len(top_samples)} SAMPLES WITH MOST COMPLETE DATA")
    print("="*80)
    
    for i, (idx, row) in enumerate(top_samples.iterrows(), 1):
        print(f"\n{i}. SRR ID: {row['srr_id']}")
        print(f"   Completeness: {int(row['completeness'])}/{len(sample_cols)} ({row['completeness_pct']:.1f}%)")
        
        # Key attributes
        key_attrs = [
            ('Sex', 'sample_sex'),
            ('Age', 'sample_age'),
            ('Disease', 'sample_disease'),
            ('Study Disease', 'sample_study_disease'),
            ('Tissue', 'sample_tissue'),
            ('Body Site', 'sample_body_site'),
            ('Cell Type', 'sample_cell_type'),
            ('Treatment', 'sample_treatment'),
            ('Genotype', 'sample_genotype'),
        ]
        
        for label, col in key_attrs:
            if col in df.columns and pd.notna(row[col]):
                val = str(row[col])[:50] + '...' if len(str(row[col])) > 50 else str(row[col])
                print(f"   {label:15s}: {val}")
        
        print(f"   Study: {row.get('study_accession', 'N/A')}")
    
    print("\n" + "="*80)
    
    return top_samples


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Find samples with complete metadata')
    parser.add_argument('--top', type=int, default=5, help='Number of top samples')
    parser.add_argument('--min-pct', type=float, help='Minimum completeness percentage')
    parser.add_argument('--require', nargs='+', help='Required fields (e.g., sex age disease)')
    parser.add_argument('--csv', type=str, help='Path to CSV file')
    parser.add_argument('--output', type=str, help='Output CSV file')
    parser.add_argument('--show-all', action='store_true', help='Show all attributes of top sample')
    
    args = parser.parse_args()
    
    # Find CSV file
    if args.csv:
        csv_file = args.csv
    else:
        metadata_dir = Path(__file__).parent.parent / 'data' / 'metadata'
        csv_files = list(metadata_dir.glob('sra_metadata_complete_*.csv'))
        if not csv_files:
            print("No metadata files found!")
            sys.exit(1)
        csv_file = max(csv_files, key=lambda p: p.stat().st_mtime)
        print(f"Using: {csv_file.name}\n")
    
    # Find complete samples
    top_samples = find_complete_samples(
        csv_file,
        top_n=args.top,
        min_completeness=args.min_pct,
        required_fields=args.require
    )
    
    if top_samples is None:
        return
    
    # Save to file
    if args.output:
        top_samples.to_csv(args.output, index=False)
        print(f"\nSaved to: {args.output}")
    
    # Show all attributes of top sample
    if args.show_all and len(top_samples) > 0:
        print("\n" + "="*80)
        print("ALL ATTRIBUTES OF TOP SAMPLE:")
        print("="*80)
        
        best = top_samples.iloc[0]
        print(f"SRR ID: {best['srr_id']}\n")
        
        sample_cols = [c for c in top_samples.columns if c.startswith('sample_')]
        filled_attrs = []
        
        for col in sample_cols:
            if pd.notna(best[col]):
                filled_attrs.append((col.replace('sample_', ''), best[col]))
        
        for attr, val in sorted(filled_attrs):
            val_str = str(val)[:70] + '...' if len(str(val)) > 70 else str(val)
            print(f"  {attr:45s}: {val_str}")
        
        print(f"\nTotal filled attributes: {len(filled_attrs)}")
        print("="*80)


if __name__ == '__main__':
    main()
