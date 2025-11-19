"""
Конфигурация для пайплайна обработки SRA данных
"""

import os
from pathlib import Path

# Базовые директории
BASE_DIR = Path(__file__).parent.parent.absolute()
DATA_DIR = BASE_DIR / 'data'
SRA_DIR = DATA_DIR / 'sra'  # Директория для SRA файлов
FASTQ_DIR = DATA_DIR / 'fastq'
QC_DIR = DATA_DIR / 'qc'
FILTERED_DIR = DATA_DIR / 'filtered'
METADATA_DIR = DATA_DIR / 'metadata'
RESULTS_DIR = DATA_DIR / 'results'

# Создаем все необходимые директории
for dir_path in [DATA_DIR, SRA_DIR, FASTQ_DIR, QC_DIR, FILTERED_DIR, METADATA_DIR, RESULTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# Параметры FastQC
FASTQC_THREADS = 4

# Параметры для Trimmomatic
TRIMMOMATIC_PARAMS = {
    'ILLUMINACLIP': str(DATA_DIR / 'adapters' / 'TruSeq3-PE.fa:2:30:10'),
    'LEADING': 3,
    'TRAILING': 3,
    'SLIDINGWINDOW': '4:15',
    'MINLEN': 36
}

# Параметры MultiQC
MULTIQC_PARAMS = {
    'force': True,
    'quiet': False
}

# Список заболеваний
DISEASES = {
    'Aortic Valve Disease': 'aortic_valve_disease',
    'Marfan Syndrome': 'marfan_syndrome',
    'Loeys-Dietz Syndrome': 'loeys_dietz_syndrome',
    'Short QT Syndrome': 'short_qt_syndrome',
    'Familial Hypercholesterolemia': 'familial_hypercholesterolemia'
}

# Максимальное количество параллельных процессов
MAX_PARALLEL_DOWNLOADS = 4
MAX_PARALLEL_PROCESSES = 4

# Пути к входным данным
CSV_FILE = DATA_DIR / 'sorted_sra_samples_by_disease_MOROZ_file.csv'

# Настройки для отчета
REPORT_TEMPLATE = """
# Анализ данных секвенирования {disease}

## Общая информация
- Всего образцов: {total_samples}
- Успешно обработано: {processed_samples}
- Дата анализа: {date}

## Статистика качества
{quality_stats}

## Статистика фильтрации
{filtering_stats}

## Метаданные
{metadata_summary}

## Графики
{plots}
"""