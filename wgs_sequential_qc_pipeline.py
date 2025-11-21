#!/usr/bin/env python3
"""
Sequential WGS QC pipeline: download, FastQC, fastp, MultiQC, cleanup
"""
import subprocess
import os
from pathlib import Path

# --- CONFIG ---
SRR_LIST = 'data/analysis_results/wgs_public_all.txt'  # список SRR
RAW_DIR = Path('data/illumina_analysis/raw_fastq')
CLEAN_DIR = Path('data/illumina_analysis/clean_fastq')
REPORT_DIR = Path('data/illumina_analysis/qc_reports')
RAW_DIR.mkdir(parents=True, exist_ok=True)
CLEAN_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

FASTQC_THREADS = 64
FASTP_THREADS = 64

with open(SRR_LIST) as f:
    srr_ids = [line.strip() for line in f if line.strip()]

for srr in srr_ids:

    # Проверка наличия отчетов
    raw_report = REPORT_DIR / f"{srr}_multiqc_raw.html"
    clean_report = REPORT_DIR / f"{srr}_multiqc_clean.html"
    if raw_report.exists() and clean_report.exists():
        print(f"[SKIP] {srr} уже обработан (отчеты найдены)")
        continue

    # Проверка свободного места на диске
    statvfs = os.statvfs(str(RAW_DIR))
    free_gb = statvfs.f_frsize * statvfs.f_bavail / (1024 ** 3)
    if free_gb < 120:
        print(f"[SKIP] {srr}: недостаточно свободного места ({free_gb:.1f} GB < 120 GB)")
        continue

    print(f"\n=== Processing {srr} ===")
    sample_raw = RAW_DIR / srr
    sample_clean = CLEAN_DIR / srr
    sample_raw.mkdir(exist_ok=True)
    sample_clean.mkdir(exist_ok=True)

    # 1. Download SRA and convert to FASTQ
    print("  Downloading SRA and converting to FASTQ...")
    subprocess.run(['prefetch', srr, '-O', str(sample_raw), '--max-size', '100G'], check=True)
    sra_path = sample_raw / srr / f"{srr}.sra"
    if not sra_path.exists():
        raise FileNotFoundError(f"SRA file not found: {sra_path}")
    subprocess.run(['fasterq-dump', str(sra_path), '--split-files', '--outdir', str(sample_raw), '--threads', '8'], check=True)

    # 2. FastQC on raw FASTQ
    print("  Running FastQC on raw FASTQ...")
    fastq_files = list(sample_raw.glob(f"{srr}_*.fastq"))
    raw_fastqc_dir = sample_raw / 'fastqc'
    raw_fastqc_dir.mkdir(exist_ok=True)
    subprocess.run(['fastqc', '-t', str(FASTQC_THREADS), '-o', str(raw_fastqc_dir)] + [str(f) for f in fastq_files], check=True)

    # 3. MultiQC raw
    print("  Running MultiQC on raw FastQC...")
    raw_multiqc_dir = sample_raw / 'multiqc'
    raw_multiqc_dir.mkdir(exist_ok=True)
    subprocess.run(['multiqc', str(raw_fastqc_dir), '-o', str(raw_multiqc_dir)], check=True)
    # Copy report
    for f in raw_multiqc_dir.glob('multiqc_report.html'):
        f.rename(REPORT_DIR / f"{srr}_multiqc_raw.html")

    # 4. fastp filtering
    print("  Running fastp filtering...")
    clean_fastq_files = []
    for fq in fastq_files:
        out_fq = sample_clean / fq.name.replace('.fastq', '.clean.fastq')
        fastp_html = REPORT_DIR / f"{srr}_{fq.name.replace('.fastq', '')}_fastp.html"
        subprocess.run([
            'fastp',
            '-i', str(fq),
            '-o', str(out_fq),
            '-w', str(FASTP_THREADS),
            '-h', str(fastp_html)
        ], check=True)
        clean_fastq_files.append(out_fq)

    # 5. FastQC on cleaned FASTQ
    print("  Running FastQC on cleaned FASTQ...")
    clean_fastqc_dir = sample_clean / 'fastqc'
    clean_fastqc_dir.mkdir(exist_ok=True)
    subprocess.run(['fastqc', '-t', str(FASTQC_THREADS), '-o', str(clean_fastqc_dir)] + [str(f) for f in clean_fastq_files], check=True)

    # 6. MultiQC cleaned
    print("  Running MultiQC on cleaned FastQC...")
    clean_multiqc_dir = sample_clean / 'multiqc'
    clean_multiqc_dir.mkdir(exist_ok=True)
    subprocess.run(['multiqc', str(clean_fastqc_dir), '-o', str(clean_multiqc_dir)], check=True)
    # Copy report
    for f in clean_multiqc_dir.glob('multiqc_report.html'):
        f.rename(REPORT_DIR / f"{srr}_multiqc_clean.html")

    # 7. Cleanup: remove raw/clean FASTQ and SRA
    print("  Cleaning up sample files...")
    for f in fastq_files + clean_fastq_files:
        if f.exists():
            f.unlink()
    if sra_path.exists():
        sra_path.unlink()
    # Remove sample directories
    for d in [sample_raw, sample_clean]:
        for sub in d.iterdir():
            if sub.is_dir():
                for sf in sub.iterdir():
                    if sf.is_file():
                        sf.unlink()
                sub.rmdir()
        d.rmdir()
    print(f"=== {srr} done ===")

print("\nAll samples processed. QC reports saved in:", REPORT_DIR)
