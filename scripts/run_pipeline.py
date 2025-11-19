"""
Батч-процессинг пайплайн для обработки SRA данных с экономией места
==================================================================

Алгоритм:
1. Читает список образцов из CSV
2. Обрабатывает их небольшими батчами (по умолчанию 5 образцов)
3. Для каждого батча:
   - Скачивает SRA файлы
   - Конвертирует в FASTQ
   - Выполняет QC анализ
   - Извлекает статистику и метрики
   - УДАЛЯЕТ RAW FASTQ файлы
   - Сохраняет только обработанные данные
4. Генерирует итоговый отчёт

Это позволяет обработать большие датасеты без заполнения диска!
"""

import os
import sys
import json
import shutil
import logging
import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from config.config import (
    CSV_FILE, DATA_DIR, SRA_DIR, FASTQ_DIR,
    QC_DIR, FILTERED_DIR, METADATA_DIR, RESULTS_DIR, DISEASES
)

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('batch_pipeline.log')
    ]
)

class BatchProcessor:
    """Батч-процессор для обработки SRA данных"""
    
    def __init__(self, batch_size: int = 5, keep_sra: bool = False, 
                 keep_raw_fastq: bool = False):
        """
        Args:
            batch_size: Количество образцов в батче
            keep_sra: Сохранять ли SRA файлы после конвертации
            keep_raw_fastq: Сохранять ли RAW FASTQ после обработки
        """
        self.batch_size = batch_size
        self.keep_sra = keep_sra
        self.keep_raw_fastq = keep_raw_fastq
        self.verify_tools()
        
        # Директория для результатов
        self.results_dir = RESULTS_DIR
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
    def verify_tools(self) -> None:
        """Проверка наличия необходимых инструментов"""
        required_tools = ['prefetch', 'fasterq-dump', 'fastqc']
        missing = []
        
        for tool in required_tools:
            try:
                subprocess.run([tool, '--version'], 
                             capture_output=True, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                missing.append(tool)
        
        if missing:
            raise RuntimeError(f"Missing tools: {', '.join(missing)}")
        
        logging.info("All required tools are available")
    
    def download_sra(self, accession: str, output_dir: Path) -> Optional[Path]:
        """Скачивание одного SRA файла"""
        try:
            logging.info(f"Downloading {accession}")
            
            # Prefetch скачивает в ~/ncbi/public/sra/ по умолчанию
            result = subprocess.run(['prefetch', accession], 
                         check=True, capture_output=True, text=True)
            
            # Ищем файл в разных возможных местах
            possible_locations = [
                Path.home() / 'ncbi' / 'public' / 'sra' / f"{accession}.sra",
                Path.home() / 'ncbi' / 'public' / 'sra' / accession / f"{accession}.sra",
                output_dir / f"{accession}.sra",
                output_dir / accession / f"{accession}.sra"
            ]
            
            for sra_file in possible_locations:
                if sra_file.exists():
                    logging.info(f"Found SRA file: {sra_file}")
                    return sra_file
            
            logging.error(f"SRA file not found in any expected location for {accession}")
            logging.error(f"Prefetch output: {result.stdout}")
            return None
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Error downloading {accession}: {e.stderr}")
            return None
    
    def convert_to_fastq(self, sra_file: Path, output_dir: Path) -> List[Path]:
        """Конвертация SRA в FASTQ"""
        try:
            accession = sra_file.stem
            logging.info(f"Converting {accession} to FASTQ")
            
            output_dir.mkdir(parents=True, exist_ok=True)
            
            subprocess.run([
                'fasterq-dump',
                '--split-files',
                '--outdir', str(output_dir),
                '--threads', '2',
                '--progress',
                str(sra_file)
            ], check=True, capture_output=True)
            
            # Находим сгенерированные FASTQ файлы
            fastq_files = list(output_dir.glob(f"{accession}*.fastq"))
            
            if not self.keep_sra:
                sra_file.unlink()
                logging.info(f"Removed SRA file: {sra_file}")
            
            return fastq_files
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error converting {sra_file}: {e.stderr.decode()}")
            return []
    
    def run_fastqc(self, fastq_file: Path, output_dir: Path) -> Dict:
        """Запуск FastQC и извлечение метрик"""
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            
            logging.info(f"Running FastQC on {fastq_file.name}")
            
            subprocess.run([
                'fastqc',
                '--outdir', str(output_dir),
                '--threads', '1',
                '--quiet',
                str(fastq_file)
            ], check=True, capture_output=True)
            
            # Извлекаем метрики из FastQC отчёта
            metrics = self.extract_fastqc_metrics(fastq_file, output_dir)
            return metrics
            
        except subprocess.CalledProcessError as e:
            logging.error(f"FastQC error: {e.stderr.decode()}")
            return {}
    
    def extract_fastqc_metrics(self, fastq_file: Path, qc_dir: Path) -> Dict:
        """Извлечение метрик из FastQC отчёта"""
        base_name = fastq_file.stem
        fastqc_data = qc_dir / f"{base_name}_fastqc" / "fastqc_data.txt"
        
        metrics = {
            'filename': fastq_file.name,
            'total_sequences': 0,
            'sequence_length': '',
            'gc_content': 0,
            'quality_score': 0
        }
        
        if not fastqc_data.exists():
            return metrics
        
        try:
            with open(fastqc_data, 'r') as f:
                for line in f:
                    if line.startswith('Total Sequences'):
                        metrics['total_sequences'] = int(line.split('\t')[1])
                    elif line.startswith('Sequence length'):
                        metrics['sequence_length'] = line.split('\t')[1].strip()
                    elif line.startswith('%GC'):
                        metrics['gc_content'] = float(line.split('\t')[1])
        except Exception as e:
            logging.error(f"Error extracting metrics: {str(e)}")
        
        return metrics
    
    def get_file_size_mb(self, file_path: Path) -> float:
        """Получение размера файла в МБ"""
        return file_path.stat().st_size / (1024 * 1024)
    
    def process_sample(self, accession: str, metadata: Dict, 
                      temp_dir: Path) -> Dict:
        """Полная обработка одного образца"""
        result = {
            'accession': accession,
            'status': 'failed',
            'error': None,
            'metrics': {},
            'metadata': metadata
        }
        
        try:
            # Используем fasterq-dump напрямую - он сам скачает SRA если нужно
            fastq_dir = temp_dir / 'fastq' / accession
            fastq_dir.mkdir(parents=True, exist_ok=True)
            
            logging.info(f"Processing {accession} (download + convert to FASTQ)")
            
            # fasterq-dump автоматически скачает SRA и сконвертирует в FASTQ
            subprocess.run([
                'fasterq-dump',
                accession,
                '--split-files',
                '--outdir', str(fastq_dir),
                '--threads', '2',
                '--progress'
            ], check=True, capture_output=True, text=True)
            
            # Находим сгенерированные FASTQ файлы
            fastq_files = list(fastq_dir.glob(f"{accession}*.fastq"))
            
            if not fastq_files:
                result['error'] = 'No FASTQ files generated'
                return result
            
            result['fastq_count'] = len(fastq_files)
            result['fastq_size_mb'] = sum(
                self.get_file_size_mb(f) for f in fastq_files
            )
            
            # 3. Контроль качества
            qc_dir = temp_dir / 'qc' / accession
            all_metrics = []
            
            for fastq_file in fastq_files:
                metrics = self.run_fastqc(fastq_file, qc_dir)
                all_metrics.append(metrics)
            
            result['metrics'] = all_metrics
            
            # 4. Сохранение QC отчётов в постоянное хранилище
            disease = metadata.get('disease', 'unknown')
            disease_slug = DISEASES.get(disease, 'unknown')
            
            permanent_qc_dir = QC_DIR / disease_slug / accession
            if qc_dir.exists():
                shutil.copytree(qc_dir, permanent_qc_dir, dirs_exist_ok=True)
            
            # 5. Удаление RAW FASTQ файлов для экономии места
            if not self.keep_raw_fastq:
                for fastq_file in fastq_files:
                    fastq_file.unlink()
                logging.info(f"Removed RAW FASTQ files for {accession}")
            
            result['status'] = 'success'
            
        except subprocess.CalledProcessError as e:
            result['error'] = f"Command failed: {e.stderr if e.stderr else 'Unknown error'}"
            logging.error(f"Error processing {accession}: {result['error']}")
        except Exception as e:
            result['error'] = str(e)
            logging.error(f"Error processing {accession}: {str(e)}")
        
        return result
    
    def process_batch(self, batch_df: pd.DataFrame, batch_num: int) -> List[Dict]:
        """Обработка одного батча образцов"""
        logging.info(f"Processing batch {batch_num} ({len(batch_df)} samples)")
        
        # Временная директория для батча
        temp_dir = DATA_DIR / f'temp_batch_{batch_num}'
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        results = []
        
        try:
            # Обработка образцов последовательно (чтобы не заполнить диск)
            for idx, row in batch_df.iterrows():
                # run_accession может содержать несколько ID через запятую
                run_accessions = row['run_accession'].strip().split(',')
                metadata = row.to_dict()
                
                for accession in run_accessions:
                    accession = accession.strip()
                    if accession:
                        result = self.process_sample(accession, metadata, temp_dir)
                        results.append(result)
                
                # Сохраняем промежуточные результаты
                self.save_batch_results(results, batch_num)
            
        finally:
            # Очищаем временную директорию
            if temp_dir.exists():
                shutil.rmtree(temp_dir)
                logging.info(f"Cleaned up temporary directory: {temp_dir}")
        
        return results
    
    def save_batch_results(self, results: List[Dict], batch_num: int) -> None:
        """Сохранение результатов батча"""
        results_file = self.results_dir / f'batch_{batch_num}_results.json'
        
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        logging.info(f"Saved batch results to {results_file}")
    
    def generate_summary_report(self, all_results: List[Dict]) -> None:
        """Генерация итогового отчёта"""
        # Подсчёт статистики
        total = len(all_results)
        successful = sum(1 for r in all_results if r['status'] == 'success')
        failed = total - successful
        
        # Создаём DataFrame для анализа
        df_list = []
        for result in all_results:
            if result['status'] == 'success':
                row = {
                    'accession': result['accession'],
                    'disease': result['metadata'].get('disease', 'unknown'),
                    'sra_size_mb': result.get('sra_size_mb', 0),
                    'fastq_size_mb': result.get('fastq_size_mb', 0),
                    'fastq_count': result.get('fastq_count', 0),
                }
                
                # Добавляем метрики из первого FASTQ файла
                if result['metrics']:
                    metrics = result['metrics'][0]
                    row.update({
                        'total_sequences': metrics.get('total_sequences', 0),
                        'gc_content': metrics.get('gc_content', 0),
                        'sequence_length': metrics.get('sequence_length', '')
                    })
                
                df_list.append(row)
        
        summary_df = pd.DataFrame(df_list)
        summary_file = self.results_dir / 'processing_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        
        # Текстовый отчёт
        report = f"""
========================================
Batch Processing Summary Report
========================================
Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Overall Statistics:
------------------
Total samples processed: {total}
Successful: {successful}
Failed: {failed}
Success rate: {successful/total*100:.1f}%

"""
        
        if not summary_df.empty:
            report += f"""
Data Statistics:
---------------
Total SRA data downloaded: {summary_df['sra_size_mb'].sum():.2f} MB
Total FASTQ data generated: {summary_df['fastq_size_mb'].sum():.2f} MB

By Disease:
-----------
"""
            disease_stats = summary_df.groupby('disease').agg({
                'accession': 'count',
                'total_sequences': 'sum',
                'sra_size_mb': 'sum'
            }).round(2)
            
            report += disease_stats.to_string()
        
        report_file = self.results_dir / 'summary_report.txt'
        with open(report_file, 'w') as f:
            f.write(report)
        
        logging.info(f"Generated summary report: {report_file}")
        print(report)

def main(disease: Optional[str] = None, batch_size: int = 5, 
         max_samples: Optional[int] = None):
    """
    Основная функция батч-процессинга
    
    Args:
        disease: Конкретное заболевание для обработки (опционально)
        batch_size: Размер батча
        max_samples: Максимальное количество образцов (для тестирования)
    """
    logging.info("Starting batch processing pipeline")
    
    # Проверка CSV файла
    if not CSV_FILE.exists():
        logging.error(f"CSV file not found: {CSV_FILE}")
        sys.exit(1)
    
    # Загрузка данных
    try:
        df = pd.read_csv(CSV_FILE, sep=';', skiprows=1)
        logging.info(f"Loaded {len(df)} samples from CSV")
    except Exception as e:
        logging.error(f"Error loading CSV: {str(e)}")
        sys.exit(1)
    
    # Фильтрация по заболеванию
    if disease:
        df = df[df['disease'] == disease]
        logging.info(f"Filtered to {len(df)} samples for disease: {disease}")
    
    # Ограничение количества образцов
    if max_samples:
        df = df.head(max_samples)
        logging.info(f"Limited to {max_samples} samples")
    
    if df.empty:
        logging.error("No samples to process")
        sys.exit(1)
    
    # Инициализация процессора
    processor = BatchProcessor(
        batch_size=batch_size,
        keep_sra=False,  # Не сохраняем SRA файлы
        keep_raw_fastq=False  # Не сохраняем RAW FASTQ
    )
    
    # Разбиваем на батчи
    all_results = []
    num_batches = (len(df) + batch_size - 1) // batch_size
    
    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(df))
        batch_df = df.iloc[start_idx:end_idx]
        
        logging.info(f"\n{'='*50}")
        logging.info(f"BATCH {batch_num + 1}/{num_batches}")
        logging.info(f"{'='*50}\n")
        
        batch_results = processor.process_batch(batch_df, batch_num + 1)
        all_results.extend(batch_results)
    
    # Генерация итогового отчёта
    processor.generate_summary_report(all_results)
    
    logging.info("Batch processing completed!")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Batch processing pipeline for SRA data'
    )
    parser.add_argument(
        '--disease',
        type=str,
        help='Specific disease to process'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=5,
        help='Number of samples per batch (default: 5)'
    )
    parser.add_argument(
        '--max-samples',
        type=int,
        help='Maximum number of samples to process (for testing)'
    )
    
    args = parser.parse_args()
    
    main(
        disease=args.disease,
        batch_size=args.batch_size,
        max_samples=args.max_samples
    )
