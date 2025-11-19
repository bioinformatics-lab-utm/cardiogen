"""
Контроль качества и фильтрация FASTQ файлов
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Add config to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config.config import (
    FASTQ_DIR, QC_DIR, FILTERED_DIR,
    TRIMMOMATIC_PARAMS, MAX_PARALLEL_PROCESSES
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('fastq_qc.log')
    ]
)

class QualityControl:
    def __init__(self):
        self.verify_tools()
        
    def verify_tools(self) -> None:
        """Проверка наличия инструментов для QC"""
        required_tools = ['fastqc', 'multiqc', 'trimmomatic']
        for tool in required_tools:
            try:
                subprocess.run([tool, '--version'], 
                             capture_output=True, check=True)
                logging.info(f"Found {tool}")
            except (subprocess.CalledProcessError, FileNotFoundError):
                raise RuntimeError(f"{tool} not found. Please install {tool}")
    
    def run_fastqc(self, fastq_file: Path) -> bool:
        """Запуск FastQC"""
        try:
            output_dir = QC_DIR / fastq_file.parent.name
            output_dir.mkdir(parents=True, exist_ok=True)
            
            subprocess.run([
                'fastqc',
                '--outdir', str(output_dir),
                '--threads', '1',
                str(fastq_file)
            ], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"FastQC error on {fastq_file}: {e.stderr.decode()}")
            return False
    
    def run_trimmomatic(self, 
                       input_files: Tuple[Path, Path],
                       output_dir: Path) -> bool:
        """Запуск Trimmomatic для парных ридов"""
        try:
            r1, r2 = input_files
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Формируем имена выходных файлов
            base_name = r1.stem.replace('_1', '')
            paired_out_1 = output_dir / f"{base_name}_1P.fastq.gz"
            paired_out_2 = output_dir / f"{base_name}_2P.fastq.gz"
            unpaired_out_1 = output_dir / f"{base_name}_1U.fastq.gz"
            unpaired_out_2 = output_dir / f"{base_name}_2U.fastq.gz"
            
            # Формируем параметры Trimmomatic
            params = [f"{k}:{v}" if k != 'ILLUMINACLIP' else v 
                     for k, v in TRIMMOMATIC_PARAMS.items()]
            
            # Определяем путь к JAR файлу Trimmomatic
            trimmomatic_jar = '/usr/share/java/trimmomatic.jar'
            if not os.path.exists(trimmomatic_jar):
                trimmomatic_jar = os.path.expanduser('~/bin/Trimmomatic/trimmomatic.jar')
                if not os.path.exists(trimmomatic_jar):
                    raise FileNotFoundError("Trimmomatic JAR file not found")
            
            cmd = [
                'java', '-jar', trimmomatic_jar,
                'PE',
                '-threads', '1',
                str(r1), str(r2),
                str(paired_out_1), str(unpaired_out_1),
                str(paired_out_2), str(unpaired_out_2),
                *params
            ]
            
            subprocess.run(cmd, check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Trimmomatic error: {e.stderr.decode()}")
            return False

    def run_multiqc(self, input_dir: Path, output_dir: Path) -> bool:
        """Запуск MultiQC для агрегации результатов"""
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            
            subprocess.run([
                'multiqc',
                str(input_dir),
                '-o', str(output_dir),
                '-f'  # Force overwrite
            ], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"MultiQC error: {e.stderr.decode()}")
            return False

def process_fastq_files(disease_dir: Path) -> None:
    """Обработка FASTQ файлов для конкретного заболевания"""
    qc = QualityControl()
    
    # Находим все пары FASTQ файлов
    fastq_pairs = {}
    for fastq in disease_dir.glob("*_1.fastq"):
        base = fastq.stem[:-2]
        pair = fastq.parent / f"{base}_2.fastq"
        if pair.exists():
            fastq_pairs[base] = (fastq, pair)
    
    # Параллельная обработка
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL_PROCESSES) as executor:
        # Запускаем FastQC
        qc_futures = []
        for pair in fastq_pairs.values():
            for fastq in pair:
                qc_futures.append(
                    executor.submit(qc.run_fastqc, fastq)
                )
        
        # Ждем завершения FastQC
        for future in as_completed(qc_futures):
            try:
                success = future.result()
                if not success:
                    logging.error("FastQC failed for some files")
            except Exception as e:
                logging.error(f"Error in FastQC: {str(e)}")
        
        # Запускаем Trimmomatic
        trim_futures = []
        filtered_dir = FILTERED_DIR / disease_dir.name
        for pair in fastq_pairs.values():
            trim_futures.append(
                executor.submit(qc.run_trimmomatic, pair, filtered_dir)
            )
        
        # Ждем завершения Trimmomatic
        for future in as_completed(trim_futures):
            try:
                success = future.result()
                if not success:
                    logging.error("Trimmomatic failed for some files")
            except Exception as e:
                logging.error(f"Error in Trimmomatic: {str(e)}")
    
    # Запускаем MultiQC для агрегации результатов
    qc_dir = QC_DIR / disease_dir.name
    qc.run_multiqc(qc_dir, qc_dir)

def main(selected_diseases: List[str] = None):
    """Основная функция пайплайна QC"""
    logging.info("Starting QC pipeline")
    
    # Check if FASTQ directory exists and has files
    if not FASTQ_DIR.exists():
        logging.error(f"FASTQ directory {FASTQ_DIR} does not exist")
        return
        
    if not any(FASTQ_DIR.iterdir()):
        logging.error(f"No files found in FASTQ directory {FASTQ_DIR}")
        return
        
    if selected_diseases is None:
        # Обрабатываем все поддиректории в FASTQ_DIR
        selected_diseases = [d.name for d in FASTQ_DIR.iterdir() if d.is_dir()]
        
    if not selected_diseases:
        logging.error("No disease directories found to process")
    
    for disease in selected_diseases:
        disease_dir = FASTQ_DIR / disease
        if disease_dir.exists():
            logging.info(f"Processing {disease} data")
            process_fastq_files(disease_dir)
        else:
            logging.error(f"Directory for {disease} not found")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Run QC on FASTQ files')
    parser.add_argument('--diseases', nargs='+', 
                       help='List of diseases to process (optional)')
    
    args = parser.parse_args()
    
    main(args.diseases)