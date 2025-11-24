#!/usr/bin/env python3
"""
Sequential WGS QC pipeline: download, FastQC, fastp, MultiQC, cleanup
Optimized for storage efficiency with resume capability
"""
import subprocess
import shutil
from pathlib import Path

# --- CONFIG ---
SRR_LIST = 'data/analysis_results/wgs_public_all.txt'
RAW_DIR = Path('data/illumina_analysis/raw_fastq')
CLEAN_DIR = Path('data/illumina_analysis/clean_fastq')
REPORT_DIR = Path('data/illumina_analysis/qc_reports')

# Centralized FastQC report directories
FASTQC_RAW_DIR = REPORT_DIR / 'fastqc_raw'
FASTQC_CLEAN_DIR = REPORT_DIR / 'fastqc_clean'

# Create all directories
for d in [RAW_DIR, CLEAN_DIR, REPORT_DIR, FASTQC_RAW_DIR, FASTQC_CLEAN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

FASTQC_THREADS = 64
FASTP_THREADS = 64

with open(SRR_LIST) as f:
    srr_ids = [line.strip() for line in f if line.strip()]

def is_sample_processed(srr):
    """Check if sample already has FastQC reports (both raw and clean)"""
    raw_fastqc_sample_dir = FASTQC_RAW_DIR / srr
    clean_fastqc_sample_dir = FASTQC_CLEAN_DIR / srr
    
    # Проверяем наличие FastQC отчетов (хотя бы один html и один zip файл в каждой директории)
    raw_has_reports = (raw_fastqc_sample_dir.exists() and 
                       list(raw_fastqc_sample_dir.glob("*_fastqc.html")) and
                       list(raw_fastqc_sample_dir.glob("*_fastqc.zip")))
    
    clean_has_reports = (clean_fastqc_sample_dir.exists() and 
                         list(clean_fastqc_sample_dir.glob("*_fastqc.html")) and
                         list(clean_fastqc_sample_dir.glob("*_fastqc.zip")))
    
    return raw_has_reports and clean_has_reports

def find_existing_sra(srr):
    """Search for existing SRA file in common locations"""
    # Поиск в текущей рабочей директории
    sample_raw = RAW_DIR / srr
    sra_path = sample_raw / srr / f"{srr}.sra"
    if sra_path.exists():
        return sra_path
    
    # Поиск в ncbi/public/sra (стандартное расположение для prefetch)
    home = Path.home()
    ncbi_sra = home / 'ncbi' / 'public' / 'sra' / f"{srr}.sra"
    if ncbi_sra.exists():
        return ncbi_sra
    
    return None

def find_existing_fastq(srr):
    """Search for existing FASTQ files"""
    sample_raw = RAW_DIR / srr
    if sample_raw.exists():
        fastq_files = list(sample_raw.glob(f"{srr}_*.fastq"))
        if fastq_files:
            return fastq_files
    return None

def process_srr(srr):
    print(f"\n=== Processing {srr} ===")
    sample_raw = RAW_DIR / srr
    sample_clean = CLEAN_DIR / srr
    fastq_files = []
    clean_fastq_files = []
    
    # Проверка: образец уже полностью обработан?
    if is_sample_processed(srr):
        print(f"  ✓ {srr} already processed, skipping...")
        return
    
    try:
        # 1. Поиск или загрузка SRA
        existing_sra = find_existing_sra(srr)
        if existing_sra:
            print(f"  ✓ Found existing SRA: {existing_sra}")
            sra_path = existing_sra
        else:
            print("  Downloading SRA...")
            sample_raw.mkdir(parents=True, exist_ok=True)
            sra_path = sample_raw / srr / f"{srr}.sra"
            subprocess.run(['prefetch', srr, '-O', str(sample_raw), '--max-size', '100G'], check=True)
            if not sra_path.exists():
                print(f"[SKIP] {srr}: SRA file not found after prefetch")
                raise Exception("SRA not found")

        # 2. Поиск или конвертация FASTQ
        existing_fastq = find_existing_fastq(srr)
        if existing_fastq:
            print(f"  ✓ Found existing FASTQ files: {len(existing_fastq)} files")
            fastq_files = existing_fastq
        else:
            print("  Converting SRA to FASTQ...")
            sample_raw.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                'fasterq-dump', str(sra_path), 
                '--split-files', 
                '--outdir', str(sample_raw), 
                '--threads', '8'
            ], check=True)
            fastq_files = list(sample_raw.glob(f"{srr}_*.fastq"))
            
            if not fastq_files:
                print(f"[SKIP] {srr}: No FASTQ files generated")
                raise Exception("FASTQ conversion failed")
        
        # Удалить SRA сразу после конвертации (или после проверки что FASTQ есть)
        if sra_path.exists():
            print("  Removing SRA file...")
            sra_path.unlink()
        # Удалить директорию SRA
        sra_dir = sample_raw / srr
        if sra_dir.exists() and sra_dir.is_dir():
            shutil.rmtree(sra_dir, ignore_errors=True)

        # 3. FastQC на raw FASTQ - проверяем, не сделан ли уже
        raw_fastqc_sample_dir = FASTQC_RAW_DIR / srr
        raw_fastqc_done = (raw_fastqc_sample_dir.exists() and 
                          list(raw_fastqc_sample_dir.glob("*_fastqc.html")))
        
        if raw_fastqc_done:
            print("  ✓ Raw FastQC already done, skipping...")
        else:
            print("  Running FastQC on raw FASTQ...")
            raw_fastqc_sample_dir.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                'fastqc', 
                '--nogroup', 
                '-t', str(FASTQC_THREADS), 
                '-o', str(raw_fastqc_sample_dir)
            ] + [str(f) for f in fastq_files], check=True)

        # 4. fastp фильтрация
        print("  Running fastp filtering...")
        sample_clean.mkdir(parents=True, exist_ok=True)
        for fq in fastq_files:
            out_fq = sample_clean / fq.name.replace('.fastq', '.clean.fastq')
            fastp_html = REPORT_DIR / f"{srr}_{fq.stem}_fastp.html"
            fastp_json = REPORT_DIR / f"{srr}_{fq.stem}_fastp.json"
            
            # Проверяем, не сделан ли уже fastp для этого файла
            if out_fq.exists() or fastp_html.exists():
                print(f"    ✓ fastp already done for {fq.name}, skipping...")
                if out_fq.exists():
                    clean_fastq_files.append(out_fq)
                continue
            
            subprocess.run([
                'fastp',
                '-i', str(fq),
                '-o', str(out_fq),
                '-w', str(FASTP_THREADS),
                '-h', str(fastp_html),
                '-j', str(fastp_json)
            ], check=True)
            clean_fastq_files.append(out_fq)
        
        # Если clean_fastq_files пуст, значит все уже было обработано - ищем существующие
        if not clean_fastq_files:
            clean_fastq_files = list(sample_clean.glob(f"{srr}_*.clean.fastq"))
        
        # Удалить raw FASTQ сразу после fastp
        print("  Removing raw FASTQ files...")
        for fq in fastq_files:
            if fq.exists():
                fq.unlink()
        fastq_files.clear()

        # 5. FastQC на очищенных FASTQ - проверяем, не сделан ли уже
        clean_fastqc_sample_dir = FASTQC_CLEAN_DIR / srr
        clean_fastqc_done = (clean_fastqc_sample_dir.exists() and 
                            list(clean_fastqc_sample_dir.glob("*_fastqc.html")))
        
        if clean_fastqc_done:
            print("  ✓ Clean FastQC already done, skipping...")
        else:
            print("  Running FastQC on clean FASTQ...")
            clean_fastqc_sample_dir.mkdir(parents=True, exist_ok=True)
            subprocess.run([
                'fastqc', 
                '--nogroup', 
                '-t', str(FASTQC_THREADS), 
                '-o', str(clean_fastqc_sample_dir)
            ] + [str(f) for f in clean_fastq_files], check=True)

        # Удалить clean FASTQ сразу после FastQC
        print("  Removing clean FASTQ files...")
        for fq in clean_fastq_files:
            if fq.exists():
                fq.unlink()
        clean_fastq_files.clear()

        print(f"  ✓ {srr} successfully processed")

    except Exception as e:
        print(f"[ERROR] {srr}: {e}")
        # В случае ошибки все равно пытаемся очистить
        for f in fastq_files + clean_fastq_files:
            if f.exists():
                try:
                    f.unlink()
                except Exception:
                    pass
    
    finally:
        # Финальная очистка: удаляем временные файлы и пустые директории
        print("  Final cleanup...")
        
        # Удалить оставшиеся .tmp и .sra файлы
        for pattern in ["*.tmp", "*.sra"]:
            if sample_raw.exists():
                for tmp in sample_raw.glob(pattern):
                    try:
                        tmp.unlink()
                    except Exception:
                        pass
        
        # Удалить директории образцов
        for d in [sample_raw, sample_clean]:
            if d.exists():
                try:
                    shutil.rmtree(d, ignore_errors=True)
                except Exception:
                    pass
        
        print(f"=== {srr} done ===\n")

# Фильтрация: обработать только необработанные образцы
print("Checking for already processed samples...")
samples_to_process = []
samples_skipped = []

for srr in srr_ids:
    if is_sample_processed(srr):
        samples_skipped.append(srr)
    else:
        samples_to_process.append(srr)

print(f"\nTotal samples: {len(srr_ids)}")
print(f"Already processed: {len(samples_skipped)}")
print(f"To process: {len(samples_to_process)}")

if samples_skipped:
    print(f"\nSkipped samples: {', '.join(samples_skipped[:10])}" + 
          (f"... and {len(samples_skipped)-10} more" if len(samples_skipped) > 10 else ""))

if not samples_to_process:
    print("\n✓ All samples already processed!")
else:
    # Запуск обработки только необработанных SRR
    print(f"\nStarting processing of {len(samples_to_process)} samples...\n")
    for i, srr in enumerate(samples_to_process, 1):
        print(f"[{i}/{len(samples_to_process)}]")
        process_srr(srr)

# Глобальные MultiQC отчеты (всегда генерируем заново)
print("\n" + "="*60)
print("Generating global MultiQC reports...")
print("="*60)

# Global MultiQC для raw данных
print("\n1. Creating global MultiQC for raw FastQC reports...")
if list(FASTQC_RAW_DIR.glob("*//*.html")):
    global_raw_multiqc = REPORT_DIR / 'multiqc_raw_global'
    global_raw_multiqc.mkdir(exist_ok=True)
    subprocess.run([
        'multiqc', 
        str(FASTQC_RAW_DIR), 
        '-o', str(global_raw_multiqc),
        '--force'
    ], check=True)
    for f in global_raw_multiqc.glob('multiqc_report.html'):
        f.rename(REPORT_DIR / "multiqc_raw_global.html")
    # Удалить временную директорию
    shutil.rmtree(global_raw_multiqc, ignore_errors=True)
    print("  ✓ Raw global MultiQC created")
else:
    print("  ⚠ No raw FastQC reports found, skipping raw global MultiQC")

# Global MultiQC для clean данных
print("\n2. Creating global MultiQC for clean FastQC reports...")
if list(FASTQC_CLEAN_DIR.glob("*//*.html")):
    global_clean_multiqc = REPORT_DIR / 'multiqc_clean_global'
    global_clean_multiqc.mkdir(exist_ok=True)
    subprocess.run([
        'multiqc', 
        str(FASTQC_CLEAN_DIR), 
        '-o', str(global_clean_multiqc),
        '--force'
    ], check=True)
    for f in global_clean_multiqc.glob('multiqc_report.html'):
        f.rename(REPORT_DIR / "multiqc_clean_global.html")
    # Удалить временную директорию
    shutil.rmtree(global_clean_multiqc, ignore_errors=True)
    print("  ✓ Clean global MultiQC created")
else:
    print("  ⚠ No clean FastQC reports found, skipping clean global MultiQC")

print("\n" + "="*60)
print("Pipeline completed successfully!")
print("="*60)
print(f"\nQC reports saved in: {REPORT_DIR}")
print(f"Raw FastQC reports: {FASTQC_RAW_DIR}")
print(f"Clean FastQC reports: {FASTQC_CLEAN_DIR}")
print("\nGenerated reports:")
print("  - Individual fastp reports: {srr}_{read}_fastp.html")
print("  - Global raw MultiQC: multiqc_raw_global.html")
print("  - Global clean MultiQC: multiqc_clean_global.html")