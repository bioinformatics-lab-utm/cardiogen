#!/usr/bin/env python3
"""
Download all possible metadata for SRR IDs from CSV files.
Uses NCBI E-utilities API to fetch comprehensive metadata.
"""

import os
import sys
import json
import time
import argparse
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Set
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config.config import BASE_DIR, DATA_DIR


class MetadataDownloader:
    """Download comprehensive SRA metadata using NCBI E-utilities."""
    
    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    ESEARCH_URL = f"{EUTILS_BASE}/esearch.fcgi"
    ESUMMARY_URL = f"{EUTILS_BASE}/esummary.fcgi"
    EFETCH_URL = f"{EUTILS_BASE}/efetch.fcgi"
    
    def __init__(self, email: str = "user@example.com", output_dir: str = None):
        """
        Initialize metadata downloader.
        
        Args:
            email: Email for NCBI API (required for polite usage)
            output_dir: Directory to save metadata files
        """
        self.email = email
        self.output_dir = Path(output_dir) if output_dir else DATA_DIR / "metadata"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Rate limiting (NCBI recommends max 3 req/sec without API key, 10 req/sec with API key)
        self.last_request_time = 0
        self.min_request_interval = 0.5  # 2 requests per second to be safe
        
    def _rate_limit(self):
        """Ensure we don't exceed NCBI rate limits."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_request_interval:
            time.sleep(self.min_request_interval - elapsed)
        self.last_request_time = time.time()
    
    def search_srr(self, srr_id: str) -> str:
        """
        Search for SRR ID and get UID.
        
        Args:
            srr_id: SRR accession number
            
        Returns:
            UID string or None
        """
        self._rate_limit()
        
        params = {
            'db': 'sra',
            'term': srr_id,
            'retmode': 'json',
            'email': self.email
        }
        
        try:
            response = requests.get(self.ESEARCH_URL, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if data.get('esearchresult', {}).get('idlist'):
                return data['esearchresult']['idlist'][0]
            return None
            
        except Exception as e:
            print(f"Error searching {srr_id}: {e}")
            return None
    
    def fetch_summary(self, uid: str) -> Dict:
        """
        Fetch summary metadata using eSummary.
        
        Args:
            uid: NCBI UID
            
        Returns:
            Dictionary with summary metadata
        """
        self._rate_limit()
        
        params = {
            'db': 'sra',
            'id': uid,
            'retmode': 'json',
            'email': self.email
        }
        
        try:
            response = requests.get(self.ESUMMARY_URL, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if 'result' in data and uid in data['result']:
                return data['result'][uid]
            return {}
            
        except Exception as e:
            print(f"Error fetching summary for UID {uid}: {e}")
            return {}
    
    def fetch_full_xml(self, srr_id: str) -> Dict:
        """
        Fetch full XML metadata using eFetch.
        
        Args:
            srr_id: SRR accession number
            
        Returns:
            Dictionary with parsed XML metadata
        """
        self._rate_limit()
        
        params = {
            'db': 'sra',
            'id': srr_id,
            'rettype': 'full',
            'retmode': 'xml',
            'email': self.email
        }
        
        try:
            response = requests.get(self.EFETCH_URL, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML
            root = ET.fromstring(response.content)
            metadata = self._parse_sra_xml(root)
            return metadata
            
        except Exception as e:
            print(f"Error fetching full XML for {srr_id}: {e}")
            return {}
    
    def _parse_sra_xml(self, root: ET.Element) -> Dict:
        """
        Parse SRA XML into structured dictionary.
        
        Args:
            root: XML root element
            
        Returns:
            Dictionary with all parsed metadata
        """
        metadata = {}
        
        # Parse EXPERIMENT_PACKAGE
        for exp_pkg in root.findall('.//EXPERIMENT_PACKAGE'):
            
            # RUN information
            run = exp_pkg.find('.//RUN')
            if run is not None:
                metadata['run'] = {
                    'accession': run.get('accession'),
                    'total_spots': run.get('total_spots'),
                    'total_bases': run.get('total_bases'),
                    'size': run.get('size'),
                    'published': run.get('published'),
                    'cluster_name': run.get('cluster_name'),
                }
                
                # Run attributes
                run_attrs = {}
                for attr in run.findall('.//RUN_ATTRIBUTE'):
                    tag = attr.find('TAG')
                    value = attr.find('VALUE')
                    if tag is not None and value is not None:
                        run_attrs[tag.text] = value.text
                metadata['run']['attributes'] = run_attrs
            
            # EXPERIMENT information
            experiment = exp_pkg.find('.//EXPERIMENT')
            if experiment is not None:
                metadata['experiment'] = {
                    'accession': experiment.get('accession'),
                    'title': self._get_text(experiment, './/TITLE'),
                    'design_description': self._get_text(experiment, './/DESIGN_DESCRIPTION'),
                    'library_name': self._get_text(experiment, './/LIBRARY_NAME'),
                    'library_strategy': self._get_text(experiment, './/LIBRARY_STRATEGY'),
                    'library_source': self._get_text(experiment, './/LIBRARY_SOURCE'),
                    'library_selection': self._get_text(experiment, './/LIBRARY_SELECTION'),
                    'library_layout': self._get_library_layout(experiment),
                    'platform': self._get_platform(experiment),
                }
                
                # Experiment attributes
                exp_attrs = {}
                for attr in experiment.findall('.//EXPERIMENT_ATTRIBUTE'):
                    tag = attr.find('TAG')
                    value = attr.find('VALUE')
                    if tag is not None and value is not None:
                        exp_attrs[tag.text] = value.text
                metadata['experiment']['attributes'] = exp_attrs
            
            # SAMPLE information
            sample = exp_pkg.find('.//SAMPLE')
            if sample is not None:
                metadata['sample'] = {
                    'accession': sample.get('accession'),
                    'title': self._get_text(sample, './/TITLE'),
                    'taxon_id': self._get_text(sample, './/TAXON_ID'),
                    'scientific_name': self._get_text(sample, './/SCIENTIFIC_NAME'),
                    'common_name': self._get_text(sample, './/COMMON_NAME'),
                    'description': self._get_text(sample, './/DESCRIPTION'),
                }
                
                # Sample attributes
                sample_attrs = {}
                for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
                    tag = attr.find('TAG')
                    value = attr.find('VALUE')
                    if tag is not None and value is not None:
                        sample_attrs[tag.text] = value.text
                metadata['sample']['attributes'] = sample_attrs
            
            # STUDY information
            study = exp_pkg.find('.//STUDY')
            if study is not None:
                metadata['study'] = {
                    'accession': study.get('accession'),
                    'title': self._get_text(study, './/STUDY_TITLE'),
                    'abstract': self._get_text(study, './/STUDY_ABSTRACT'),
                    'description': self._get_text(study, './/STUDY_DESCRIPTION'),
                    'study_type': self._get_text(study, './/STUDY_TYPE'),
                }
                
                # Study attributes
                study_attrs = {}
                for attr in study.findall('.//STUDY_ATTRIBUTE'):
                    tag = attr.find('TAG')
                    value = attr.find('VALUE')
                    if tag is not None and value is not None:
                        study_attrs[tag.text] = value.text
                metadata['study']['attributes'] = study_attrs
            
            # SUBMISSION information
            submission = exp_pkg.find('.//SUBMISSION')
            if submission is not None:
                metadata['submission'] = {
                    'accession': submission.get('accession'),
                    'submission_date': submission.get('submission_date'),
                    'lab_name': submission.get('lab_name'),
                    'center_name': submission.get('center_name'),
                }
            
            # Organization
            organization = exp_pkg.find('.//Organization')
            if organization is not None:
                metadata['organization'] = {
                    'type': organization.get('type'),
                    'name': self._get_text(organization, './/Name'),
                }
        
        return metadata
    
    def _get_text(self, element: ET.Element, path: str) -> str:
        """Safely get text from XML element."""
        found = element.find(path)
        return found.text if found is not None and found.text else ""
    
    def _get_library_layout(self, experiment: ET.Element) -> str:
        """Get library layout (SINGLE or PAIRED)."""
        if experiment.find('.//PAIRED') is not None:
            return 'PAIRED'
        elif experiment.find('.//SINGLE') is not None:
            return 'SINGLE'
        return ''
    
    def _get_platform(self, experiment: ET.Element) -> Dict:
        """Get platform information."""
        platform_info = {}
        platform = experiment.find('.//PLATFORM')
        if platform is not None:
            for child in platform:
                platform_info['type'] = child.tag
                platform_info['instrument_model'] = self._get_text(child, './/INSTRUMENT_MODEL')
        return platform_info
    
    def download_metadata(self, srr_id: str) -> Dict:
        """
        Download complete metadata for a single SRR ID.
        
        Args:
            srr_id: SRR accession number
            
        Returns:
            Dictionary with all metadata
        """
        print(f"Downloading metadata for {srr_id}...")
        
        metadata = {
            'srr_id': srr_id,
            'download_timestamp': datetime.now().isoformat(),
            'summary': {},
            'full': {}
        }
        
        # Get UID
        uid = self.search_srr(srr_id)
        if uid:
            metadata['uid'] = uid
            # Get summary
            metadata['summary'] = self.fetch_summary(uid)
        
        # Get full XML metadata
        metadata['full'] = self.fetch_full_xml(srr_id)
        
        return metadata
    
    def download_batch(self, srr_ids: List[str], max_workers: int = 2, 
                      output_dir: Path = None) -> List[Dict]:
        """
        Download metadata for multiple SRR IDs in parallel.
        Saves after each record to avoid data loss.
        
        Args:
            srr_ids: List of SRR accession numbers
            max_workers: Maximum number of parallel downloads (default: 2 to avoid rate limits)
            output_dir: Directory to save files
            
        Returns:
            List of metadata dictionaries
        """
        all_metadata = []
        temp_file = output_dir / "sra_metadata_PARTIAL_IN_PROGRESS.json"
        
        # Load existing partial data if available
        if temp_file.exists():
            try:
                with open(temp_file, 'r') as f:
                    all_metadata = json.load(f)
                print(f"Resuming from {len(all_metadata)} existing records")
            except:
                pass
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_srr = {
                executor.submit(self.download_metadata, srr_id): srr_id 
                for srr_id in srr_ids
            }
            
            for i, future in enumerate(as_completed(future_to_srr), 1):
                srr_id = future_to_srr[future]
                try:
                    metadata = future.result()
                    all_metadata.append(metadata)
                    print(f"Progress: {i}/{len(srr_ids)} - {srr_id} completed")
                    
                    # Save after EVERY record
                    with open(temp_file, 'w') as f:
                        json.dump(all_metadata, f, indent=2)
                        
                except Exception as e:
                    print(f"Error downloading {srr_id}: {e}")
        
        return all_metadata


def extract_srr_ids_from_csv(csv_file: Path) -> Set[str]:
    """
    Extract all unique SRR IDs from a CSV file.
    
    Args:
        csv_file: Path to CSV file
        
    Returns:
        Set of unique SRR IDs
    """
    srr_ids = set()
    
    try:
        # Try reading with different parameters
        if 'MOROZ' in csv_file.name:
            df = pd.read_csv(csv_file, sep=';', skiprows=1)
        else:
            df = pd.read_csv(csv_file)
        
        # Find column with run accessions
        run_col = None
        for col in df.columns:
            if 'run' in col.lower() or 'srr' in col.lower() or 'err' in col.lower():
                run_col = col
                break
        
        if run_col:
            for value in df[run_col].dropna():
                # Handle comma-separated values
                if ',' in str(value):
                    for srr in str(value).split(','):
                        srr = srr.strip()
                        if srr.startswith(('SRR', 'ERR', 'DRR')):
                            srr_ids.add(srr)
                else:
                    srr = str(value).strip()
                    if srr.startswith(('SRR', 'ERR', 'DRR')):
                        srr_ids.add(srr)
        
        print(f"Extracted {len(srr_ids)} unique SRR IDs from {csv_file.name}")
        
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
    
    return srr_ids


def save_metadata(metadata_list: List[Dict], output_dir: Path):
    """
    Save metadata to JSON and CSV files.
    
    Args:
        metadata_list: List of metadata dictionaries
        output_dir: Directory to save files
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save full metadata as JSON
    json_file = output_dir / f"sra_metadata_full_{timestamp}.json"
    with open(json_file, 'w') as f:
        json.dump(metadata_list, f, indent=2)
    print(f"Saved full metadata to {json_file}")
    
    # Create flattened CSV for easy analysis
    csv_data = []
    for meta in metadata_list:
        row = {
            'srr_id': meta['srr_id'],
            'download_timestamp': meta['download_timestamp'],
        }
        
        # Add summary fields
        if 'summary' in meta and 'expxml' in meta['summary']:
            expxml = meta['summary']['expxml']
            row['summary_runs'] = meta['summary'].get('runs', '')
            row['summary_total_spots'] = meta['summary'].get('total_spots', '')
            row['summary_total_size'] = meta['summary'].get('total_size', '')
        
        # Add full metadata fields
        if 'full' in meta:
            full = meta['full']
            
            # Run info
            if 'run' in full:
                row['run_accession'] = full['run'].get('accession', '')
                row['run_total_spots'] = full['run'].get('total_spots', '')
                row['run_total_bases'] = full['run'].get('total_bases', '')
                row['run_size'] = full['run'].get('size', '')
                row['run_published'] = full['run'].get('published', '')
                
                # Run attributes (quality scores, etc.)
                if 'attributes' in full['run']:
                    for key, value in full['run']['attributes'].items():
                        row[f'run_attr_{key}'] = value
            
            # Experiment info
            if 'experiment' in full:
                exp = full['experiment']
                row['experiment_accession'] = exp.get('accession', '')
                row['experiment_title'] = exp.get('title', '')
                row['library_strategy'] = exp.get('library_strategy', '')
                row['library_source'] = exp.get('library_source', '')
                row['library_selection'] = exp.get('library_selection', '')
                row['library_layout'] = exp.get('library_layout', '')
                row['library_name'] = exp.get('library_name', '')
                row['design_description'] = exp.get('design_description', '')
                
                if 'platform' in exp:
                    row['platform_type'] = exp['platform'].get('type', '')
                    row['instrument_model'] = exp['platform'].get('instrument_model', '')
                
                # Experiment attributes
                if 'attributes' in exp:
                    for key, value in exp['attributes'].items():
                        row[f'exp_attr_{key}'] = value
            
            # Sample info - THIS IS WHERE PATIENT DATA IS!
            if 'sample' in full:
                sample = full['sample']
                row['sample_accession'] = sample.get('accession', '')
                row['sample_title'] = sample.get('title', '')
                row['organism'] = sample.get('scientific_name', '')
                row['taxon_id'] = sample.get('taxon_id', '')
                row['sample_description'] = sample.get('description', '')
                
                # Sample attributes - PATIENT METADATA (age, sex, disease, tissue, etc.)
                if 'attributes' in sample:
                    for key, value in sample['attributes'].items():
                        # Clean up key name for CSV column
                        clean_key = key.lower().replace(' ', '_').replace('-', '_')
                        row[f'sample_{clean_key}'] = value
            
            # Study info
            if 'study' in full:
                study = full['study']
                row['study_accession'] = study.get('accession', '')
                row['study_title'] = study.get('title', '')
                row['study_abstract'] = study.get('abstract', '')
                row['study_description'] = study.get('description', '')
                row['study_type'] = study.get('study_type', '')
                
                # Study attributes
                if 'attributes' in study:
                    for key, value in study['attributes'].items():
                        row[f'study_attr_{key}'] = value
            
            # Submission info
            if 'submission' in full:
                sub = full['submission']
                row['submission_accession'] = sub.get('accession', '')
                row['submission_date'] = sub.get('submission_date', '')
                row['center_name'] = sub.get('center_name', '')
                row['lab_name'] = sub.get('lab_name', '')
            
            # Organization info
            if 'organization' in full:
                org = full['organization']
                row['organization_type'] = org.get('type', '')
                row['organization_name'] = org.get('name', '')
        
        csv_data.append(row)
    
    # Save to CSV
    csv_file = output_dir / f"sra_metadata_complete_{timestamp}.csv"
    df = pd.DataFrame(csv_data)
    df.to_csv(csv_file, index=False)
    
    # Print detailed summary
    print(f"\n{'='*70}")
    print("SAVED: Complete metadata CSV")
    print(f"{'='*70}")
    print(f"File: {csv_file}")
    print(f"Size: {csv_file.stat().st_size / 1024 / 1024:.2f} MB")
    print(f"\nTotal rows: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    
    # Show column categories
    print(f"\n{'─'*70}")
    print("COLUMN CATEGORIES:")
    print(f"{'─'*70}")
    
    basic_cols = [c for c in df.columns if c in ['srr_id', 'download_timestamp']]
    run_cols = [c for c in df.columns if c.startswith('run_')]
    exp_cols = [c for c in df.columns if c.startswith('experiment_') or c.startswith('exp_attr_')]
    sample_cols = [c for c in df.columns if c.startswith('sample_')]
    study_cols = [c for c in df.columns if c.startswith('study_')]
    other_cols = [c for c in df.columns if c.startswith(('submission_', 'platform_', 'instrument_', 'library_', 'organism', 'center_', 'organization_'))]
    
    print(f"  Basic info:        {len(basic_cols)} columns")
    print(f"  Run data:          {len(run_cols)} columns (sequencing metrics)")
    print(f"  Experiment data:   {len(exp_cols)} columns (library info)")
    print(f"  Sample data:       {len(sample_cols)} columns ⭐ PATIENT/CLINICAL DATA")
    print(f"  Study data:        {len(study_cols)} columns")
    print(f"  Other metadata:    {len(other_cols)} columns")
    
    # Show patient/clinical columns
    if sample_cols:
        print(f"\n{'─'*70}")
        print("PATIENT/CLINICAL DATA COLUMNS (sample_*):")
        print(f"{'─'*70}")
        
        # Group by common keywords
        clinical_keywords = {
            'Demographics': ['sex', 'gender', 'age'],
            'Disease/Phenotype': ['disease', 'phenotype', 'affected', 'diagnosis', 'condition'],
            'Tissue/Cell': ['tissue', 'body_site', 'cell_type', 'cell_line', 'source'],
            'Genetics': ['genotype', 'mutation', 'variant', 'allele'],
            'Treatment': ['treatment', 'therapy', 'drug', 'dose'],
            'Clinical': ['clinical', 'measurement', 'lab', 'test', 'score'],
        }
        
        for category, keywords in clinical_keywords.items():
            matching = [c for c in sample_cols if any(kw in c.lower() for kw in keywords)]
            if matching:
                print(f"\n  {category}:")
                for col in matching[:5]:  # Show first 5
                    # Show example value
                    example = df[col].dropna().iloc[0] if len(df[col].dropna()) > 0 else "N/A"
                    example_str = str(example)[:40] + "..." if len(str(example)) > 40 else str(example)
                    print(f"    • {col.replace('sample_', '')}: {example_str}")
                if len(matching) > 5:
                    print(f"    ... and {len(matching) - 5} more")
        
        # Show other sample columns
        categorized = []
        for keywords in clinical_keywords.values():
            categorized.extend([c for c in sample_cols if any(kw in c.lower() for kw in keywords)])
        uncategorized = [c for c in sample_cols if c not in categorized]
        
        if uncategorized:
            print(f"\n  Other attributes: {len(uncategorized)} columns")
            for col in uncategorized[:3]:
                example = df[col].dropna().iloc[0] if len(df[col].dropna()) > 0 else "N/A"
                example_str = str(example)[:40] + "..." if len(str(example)) > 40 else str(example)
                print(f"    • {col.replace('sample_', '')}: {example_str}")
            if len(uncategorized) > 3:
                print(f"    ... and {len(uncategorized) - 3} more")
    
    # Show data completeness
    print(f"\n{'─'*70}")
    print("DATA COMPLETENESS:")
    print(f"{'─'*70}")
    
    if sample_cols:
        # Calculate % of non-null values for key patient columns
        key_patient_cols = [c for c in sample_cols if any(kw in c.lower() for kw in ['sex', 'age', 'disease', 'tissue', 'body_site'])]
        for col in key_patient_cols[:10]:
            completeness = (df[col].notna().sum() / len(df)) * 100
            print(f"  {col.replace('sample_', '')}: {completeness:.1f}% ({df[col].notna().sum()}/{len(df)} samples)")
    
    print(f"\n{'='*70}")
    
    return json_file, csv_file


def main():
    parser = argparse.ArgumentParser(
        description='Download all possible metadata for SRR IDs from CSV files'
    )
    parser.add_argument(
        '--email',
        type=str,
        default='user@example.com',
        help='Email for NCBI API (recommended for polite usage)'
    )
    parser.add_argument(
        '--csv-files',
        type=str,
        nargs='+',
        help='Specific CSV files to process (default: all CSV files in data/)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        help='Output directory for metadata files'
    )
    parser.add_argument(
        '--max-workers',
        type=int,
        default=2,
        help='Maximum number of parallel downloads (default: 2, reduce if hitting rate limits)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='Limit number of SRR IDs to process (for testing)'
    )
    
    args = parser.parse_args()
    
    # Find CSV files
    if args.csv_files:
        csv_files = [Path(f) for f in args.csv_files]
    else:
        csv_files = list(DATA_DIR.glob('*.csv'))
    
    if not csv_files:
        print("No CSV files found!")
        return
    
    print(f"Found {len(csv_files)} CSV files:")
    for csv_file in csv_files:
        print(f"  - {csv_file.name}")
    
    # Extract all unique SRR IDs
    print("\nExtracting SRR IDs from CSV files...")
    all_srr_ids = set()
    for csv_file in csv_files:
        srr_ids = extract_srr_ids_from_csv(csv_file)
        all_srr_ids.update(srr_ids)
    
    print(f"\nTotal unique SRR IDs: {len(all_srr_ids)}")
    
    # Apply limit if specified
    if args.limit:
        all_srr_ids = list(all_srr_ids)[:args.limit]
        print(f"Limited to {len(all_srr_ids)} SRR IDs for testing")
    
    # Initialize downloader
    downloader = MetadataDownloader(
        email=args.email,
        output_dir=args.output_dir
    )
    
    # Download metadata
    print(f"\nStarting metadata download with {args.max_workers} parallel workers...")
    print(f"This will take approximately {len(all_srr_ids) * 0.5 / 60:.1f} minutes")
    
    start_time = time.time()
    metadata_list = downloader.download_batch(
        list(all_srr_ids),
        max_workers=args.max_workers,
        output_dir=downloader.output_dir
    )
    elapsed = time.time() - start_time
    
    # Remove partial file after successful completion
    temp_file = downloader.output_dir / "sra_metadata_PARTIAL_IN_PROGRESS.json"
    if temp_file.exists():
        temp_file.unlink()
    
    print(f"\nDownload completed in {elapsed/60:.1f} minutes")
    print(f"Successfully downloaded metadata for {len(metadata_list)} SRR IDs")
    
    # Save metadata
    print("\n" + "="*70)
    print("SAVING RESULTS...")
    print("="*70)
    json_file, csv_file = save_metadata(metadata_list, downloader.output_dir)
    
    print("\n" + "="*70)
    print("✓ METADATA DOWNLOAD COMPLETE!")
    print("="*70)
    print(f"Processed: {len(metadata_list)} SRR IDs")
    print(f"Time elapsed: {elapsed/60:.1f} minutes")
    print(f"\nOUTPUT FILES:")
    print(f"  1. Full metadata (JSON): {json_file}")
    print(f"  2. Complete table (CSV): {csv_file} ⭐ MAIN FILE")
    print(f"\nNEXT STEPS:")
    print(f"  • Open CSV in Excel/LibreOffice to view all data")
    print(f"  • Use pandas to analyze: pd.read_csv('{csv_file.name}')")
    print(f"  • Filter by patient attributes (sex, age, disease, etc.)")
    print("="*70)


if __name__ == "__main__":
    main()
