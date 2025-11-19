#!/usr/bin/env python3
"""
View current partial metadata that's being downloaded.
"""

import json
import pandas as pd
from pathlib import Path

partial_file = Path("/home/nicolaedrabcinski/cardiogen/data/metadata/sra_metadata_PARTIAL_IN_PROGRESS.json")

if not partial_file.exists():
    print("No partial metadata file found!")
    print("Either download hasn't started or already completed.")
    exit(1)

print("Loading partial metadata...")
with open(partial_file) as f:
    data = json.load(f)

print(f"\nCurrent progress: {len(data)} records downloaded")

# Convert to simple stats
print("\n" + "="*60)
print("CURRENT DATA OVERVIEW")
print("="*60)

# Count records with patient data
patient_data_count = 0
for record in data:
    if 'full' in record and 'sample' in record['full']:
        if 'attributes' in record['full']['sample'] and record['full']['sample']['attributes']:
            patient_data_count += 1

print(f"Records with patient attributes: {patient_data_count}/{len(data)}")

# Extract some sample info
print("\n" + "="*60)
print("SAMPLE DATA")
print("="*60)

for i, record in enumerate(data[:3]):
    print(f"\nSample {i+1}: {record['srr_id']}")
    if 'full' in record and 'sample' in record['full']:
        sample = record['full']['sample']
        if 'attributes' in sample:
            attrs = sample['attributes']
            for key, val in list(attrs.items())[:5]:
                print(f"  {key}: {str(val)[:50]}")

print("\n" + "="*60)
print(f"\nTotal records: {len(data)}")
print("File will be converted to CSV when download completes.")
print("="*60)
