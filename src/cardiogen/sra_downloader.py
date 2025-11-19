"""
Загрузка и подготовка FASTQ файлов из SRA
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import time

import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add config to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from config.config import (
    FASTQ_DIR, METADATA_DIR, DISEASES, SRA_DIR,
    MAX_PARALLEL_DOWNLOADS
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('sra_downloader.log')
    ]
)

class SRADownloader:
    def __init__(self, max_retries: int = 3, wait_time: int = 10):
        self.max_retries = max_retries
        self.wait_time = wait_time
        self.verify_tools()
        
    def verify_tools(self) -> None:
        """Проверка наличия необходимых инструментов"""
        required_tools = ['prefetch', 'fasterq-dump', 'fastq-dump']
        for tool in required_tools:
            try:
                subprocess.run([tool, '--version'], 
                             capture_output=True, check=True)
                logging.info(f"Found {tool}")
            except (subprocess.CalledProcessError, FileNotFoundError):
                raise RuntimeError(f"{tool} not found. Please install SRA Toolkit")

    def download_single_run(self, run_id: str, metadata: Dict, 
                          output_dir: Path) -> bool:
        """Загрузка одного SRA запуска"""
        try:
            # Prefetch для загрузки sra файла
            logging.info(f"Prefetching {run_id}")
            subprocess.run(['prefetch', run_id], check=True,
                         capture_output=True)
            
            # Конвертация в FASTQ с сохранением парных ридов
            logging.info(f"Converting {run_id} to FASTQ")
            subprocess.run([
                'fasterq-dump',
                '--split-files',
                '--outdir', str(output_dir),
                '--threads', '4',
                run_id
            ], check=True, capture_output=True)
            
            # Сохранение метаданных
            self._save_metadata(run_id, metadata, output_dir)
            
            return True
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {run_id}: {e.stderr.decode()}")
            return False
            
    def _save_metadata(self, run_id: str, metadata: Dict, 
                      output_dir: Path) -> None:
        """Сохранение метаданных в JSON"""
        metadata_file = output_dir / f"{run_id}_metadata.json"
        pd.Series(metadata).to_json(metadata_file)

class FastqProcessor:
    def __init__(self):
        self.verify_tools()
        
    def verify_tools(self) -> None:
        """Проверка наличия инструментов для обработки FASTQ"""
        required_tools = ['fastqc', 'multiqc']
        for tool in required_tools:
            try:
                subprocess.run([tool, '--version'], 
                             capture_output=True, check=True)
                logging.info(f"Found {tool}")
            except (subprocess.CalledProcessError, FileNotFoundError):
                raise RuntimeError(f"{tool} not found. Please install {tool}")
    
    def run_fastqc(self, fastq_file: Path) -> bool:
        """Запуск FastQC для контроля качества"""
        try:
            subprocess.run([
                'fastqc',
                '--outdir', str(fastq_file.parent),
                '--threads', '1',
                str(fastq_file)
            ], check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"FastQC error on {fastq_file}: {e.stderr.decode()}")
            return False

def download_sra_files(accessions: List[str], output_dir: Path) -> List[Path]:
    """
    Скачивание SRA файлов
    
    Args:
        accessions: Список SRA идентификаторов
        output_dir: Директория для сохранения файлов
        
    Returns:
        Список путей к скачанным SRA файлам
    """
    downloader = SRADownloader()
    downloaded_files = []
    
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL_DOWNLOADS) as executor:
        future_to_acc = {
            executor.submit(
                downloader.download_single_run,
                acc,
                {'run_accession': acc},
                output_dir
            ): acc for acc in accessions
        }
        
        for future in as_completed(future_to_acc):
            acc = future_to_acc[future]
            try:
                success = future.result()
                if success:
                    sra_file = SRA_DIR / acc / f"{acc}.sra"
                    downloaded_files.append(sra_file)
                    logging.info(f"Successfully downloaded {acc}")
                else:
                    logging.error(f"Failed to download {acc}")
            except Exception as e:
                logging.error(f"Error downloading {acc}: {str(e)}")
    
    return downloaded_files

def convert_to_fastq(sra_files: List[Path], output_dir: Path) -> List[Path]:
    """
    Конвертация SRA файлов в FASTQ формат
    
    Args:
        sra_files: Список путей к SRA файлам
        output_dir: Директория для сохранения FASTQ файлов
        
    Returns:
        Список путей к сгенерированным FASTQ файлам
    """
    fastq_files = []
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for sra_file in sra_files:
        try:
            # Получаем SRA идентификатор из имени файла
            sra_id = sra_file.stem
            
            # Запускаем fasterq-dump
            logging.info(f"Converting {sra_id} to FASTQ")
            subprocess.run([
                'fasterq-dump',
                '--split-files',
                '--outdir', str(output_dir),
                '--threads', '4',
                str(sra_file)
            ], check=True, capture_output=True)
            
            # Добавляем сгенерированные файлы в список
            for fastq in output_dir.glob(f"{sra_id}*.fastq"):
                fastq_files.append(fastq)
                logging.info(f"Generated {fastq.name}")
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Error converting {sra_file}: {e.stderr.decode()}")
            continue
        except Exception as e:
            logging.error(f"Unexpected error converting {sra_file}: {str(e)}")
            continue
    
    return fastq_files

def process_disease_data(csv_file: Path, disease: str) -> None:
    """Обработка данных для конкретного заболевания"""
    # Чтение CSV и фильтрация по заболеванию
    df = pd.read_csv(csv_file, sep=';', skiprows=1)
    disease_df = df[df['disease'] == disease]
    
    # Создание директорий
    disease_dir = FASTQ_DIR / DISEASES[disease]
    disease_dir.mkdir(parents=True, exist_ok=True)
    
    metadata_dir = METADATA_DIR / DISEASES[disease]
    metadata_dir.mkdir(parents=True, exist_ok=True)
    
    # Инициализация загрузчика и процессора
    downloader = SRADownloader()
    processor = FastqProcessor()
    
    # Обработка каждого образца
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL_DOWNLOADS) as executor:
        futures = []
        for _, row in disease_df.iterrows():
            run_ids = [rid.strip() for rid in row['run_accession'].split(',')]
            for run_id in run_ids:
                metadata = row.to_dict()
                futures.append(
                    executor.submit(
                        downloader.download_single_run,
                        run_id,
                        metadata,
                        disease_dir
                    )
                )
        
        # Сбор результатов
        for future in as_completed(futures):
            try:
                success = future.result()
                if success:
                    logging.info("Successfully processed run")
                else:
                    logging.error("Failed to process run")
            except Exception as e:
                logging.error(f"Error in processing: {str(e)}")

def main(csv_file: Path, selected_diseases: Optional[List[str]] = None):
    """Основная функция пайплайна"""
    if selected_diseases is None:
        selected_diseases = list(DISEASES.keys())
    
    for disease in selected_diseases:
        logging.info(f"Processing {disease} data")
        process_disease_data(csv_file, disease)
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Download and process SRA data')
    parser.add_argument('--diseases', nargs='+', 
                       help='List of diseases to process (optional)')
    parser.add_argument('--csv', type=str, default='sorted_sra_samples_by_disease_MOROZ_file.csv',
                       help='Path to input CSV file')
    
    args = parser.parse_args()
    
    main(Path(args.csv), args.diseases)