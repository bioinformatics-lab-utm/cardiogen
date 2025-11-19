#!/usr/bin/env python3
"""
Детальный EDA для Illumina образцов
Фильтрует только Illumina платформу и проводит углубленный анализ метаданных
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Настройка стиля графиков
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 10


class IlluminaEDA:
    """Детальный EDA для Illumina образцов"""
    
    def __init__(self, metadata_file: str):
        """
        Args:
            metadata_file: Путь к CSV файлу с метаданными
        """
        print(f"Загрузка метаданных из: {metadata_file}")
        self.df_full = pd.read_csv(metadata_file, low_memory=False)
        print(f"Всего образцов в файле: {len(self.df_full)}")
        
        # Фильтруем только Illumina
        self.df = self._filter_illumina()
        print(f"Illumina образцов: {len(self.df)}")
        
        # Создаем директорию для результатов
        self.output_dir = Path('data/metadata/illumina_eda')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Определяем колонки с атрибутами пациентов
        self.sample_cols = [c for c in self.df.columns if c.startswith('sample_')]
        print(f"Атрибутов пациентов: {len(self.sample_cols)}\n")
    
    def _filter_illumina(self) -> pd.DataFrame:
        """Фильтрация только Illumina образцов"""
        # Проверяем разные поля, где может быть указана платформа
        illumina_mask = (
            (self.df_full['platform_type'].str.contains('ILLUMINA', case=False, na=False)) |
            (self.df_full['platform_type'].str.contains('Illumina', case=False, na=False)) |
            (self.df_full['instrument_model'].str.contains('Illumina', case=False, na=False))
        )
        
        return self.df_full[illumina_mask].copy()
    
    def platform_overview(self):
        """Обзор платформ Illumina"""
        print("="*80)
        print("ОБЗОР ILLUMINA ПЛАТФОРМ")
        print("="*80)
        
        # Подсчет по моделям инструментов
        instruments = self.df['instrument_model'].value_counts()
        print("\nМодели инструментов Illumina:")
        for instrument, count in instruments.items():
            pct = (count / len(self.df)) * 100
            print(f"  {instrument:40s}: {count:6d} ({pct:5.1f}%)")
        
        # График распределения инструментов
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Топ 10 инструментов
        top_instruments = instruments.head(10)
        ax1.barh(range(len(top_instruments)), top_instruments.values)
        ax1.set_yticks(range(len(top_instruments)))
        ax1.set_yticklabels(top_instruments.index)
        ax1.set_xlabel('Количество образцов')
        ax1.set_title('Топ 10 моделей Illumina инструментов')
        ax1.invert_yaxis()
        
        for i, v in enumerate(top_instruments.values):
            ax1.text(v, i, f' {v}', va='center')
        
        # Pie chart для топ 5
        top5 = instruments.head(5)
        other = instruments[5:].sum()
        if other > 0:
            plot_data = pd.concat([top5, pd.Series({'Other': other})])
        else:
            plot_data = top5
        
        ax2.pie(plot_data.values, labels=plot_data.index, autopct='%1.1f%%', startangle=90)
        ax2.set_title('Доля основных Illumina инструментов')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'illumina_instruments.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"\n✓ Сохранено: {self.output_dir / 'illumina_instruments.png'}")
    
    def library_strategy_analysis(self):
        """Анализ стратегий библиотек"""
        print("\n" + "="*80)
        print("АНАЛИЗ СТРАТЕГИЙ БИБЛИОТЕК")
        print("="*80)
        
        strategies = self.df['library_strategy'].value_counts()
        print("\nСтратегии библиотек:")
        for strategy, count in strategies.items():
            pct = (count / len(self.df)) * 100
            print(f"  {strategy:30s}: {count:6d} ({pct:5.1f}%)")
        
        # График
        fig, ax = plt.subplots(figsize=(12, 6))
        strategies.plot(kind='bar', ax=ax)
        ax.set_xlabel('Стратегия библиотеки')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение стратегий библиотек (Illumina)')
        plt.xticks(rotation=45, ha='right')
        
        # Добавляем значения на столбцы
        for i, v in enumerate(strategies.values):
            ax.text(i, v, f'{v}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'library_strategies.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'library_strategies.png'}")
    
    def library_layout_analysis(self):
        """Анализ типов лэйаута библиотек"""
        print("\n" + "="*80)
        print("АНАЛИЗ ТИПОВ ЛЭЙАУТА")
        print("="*80)
        
        layouts = self.df['library_layout'].value_counts()
        print("\nТипы лэйаута:")
        for layout, count in layouts.items():
            pct = (count / len(self.df)) * 100
            print(f"  {layout:20s}: {count:6d} ({pct:5.1f}%)")
        
        # График
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = sns.color_palette('Set2', n_colors=len(layouts))
        ax.pie(layouts.values, labels=layouts.index, autopct='%1.1f%%', 
               startangle=90, colors=colors, textprops={'fontsize': 12})
        ax.set_title('Распределение типов лэйаута библиотек (Illumina)', fontsize=14)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'library_layouts.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'library_layouts.png'}")
    
    def sequencing_depth_analysis(self):
        """Анализ глубины секвенирования"""
        print("\n" + "="*80)
        print("АНАЛИЗ ГЛУБИНЫ СЕКВЕНИРОВАНИЯ")
        print("="*80)
        
        # Конвертируем в числовой формат
        self.df['run_total_spots_num'] = pd.to_numeric(self.df['run_total_spots'], errors='coerce')
        self.df['run_total_bases_num'] = pd.to_numeric(self.df['run_total_bases'], errors='coerce')
        
        # Статистика по spots
        spots_stats = self.df['run_total_spots_num'].describe()
        print("\nСтатистика по количеству spots:")
        print(f"  Среднее:     {spots_stats['mean']:,.0f}")
        print(f"  Медиана:     {spots_stats['50%']:,.0f}")
        print(f"  Мин:         {spots_stats['min']:,.0f}")
        print(f"  Макс:        {spots_stats['max']:,.0f}")
        print(f"  Ст. откл.:   {spots_stats['std']:,.0f}")
        
        # Статистика по bases
        bases_stats = self.df['run_total_bases_num'].describe()
        print("\nСтатистика по количеству bases:")
        print(f"  Среднее:     {bases_stats['mean']:,.0f}")
        print(f"  Медиана:     {bases_stats['50%']:,.0f}")
        print(f"  Мин:         {bases_stats['min']:,.0f}")
        print(f"  Макс:        {bases_stats['max']:,.0f}")
        print(f"  Ст. откл.:   {bases_stats['std']:,.0f}")
        
        # Графики
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Histogram spots
        ax = axes[0, 0]
        self.df['run_total_spots_num'].dropna().hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Total Spots')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение Total Spots')
        ax.axvline(spots_stats['mean'], color='red', linestyle='--', label=f"Mean: {spots_stats['mean']:,.0f}")
        ax.axvline(spots_stats['50%'], color='green', linestyle='--', label=f"Median: {spots_stats['50%']:,.0f}")
        ax.legend()
        
        # Log-scale histogram spots
        ax = axes[0, 1]
        self.df['run_total_spots_num'].dropna().apply(np.log10).hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Log10(Total Spots)')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение Total Spots (log-scale)')
        
        # Histogram bases
        ax = axes[1, 0]
        self.df['run_total_bases_num'].dropna().hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Total Bases')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение Total Bases')
        ax.axvline(bases_stats['mean'], color='red', linestyle='--', label=f"Mean: {bases_stats['mean']:,.0f}")
        ax.axvline(bases_stats['50%'], color='green', linestyle='--', label=f"Median: {bases_stats['50%']:,.0f}")
        ax.legend()
        
        # Log-scale histogram bases
        ax = axes[1, 1]
        self.df['run_total_bases_num'].dropna().apply(np.log10).hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Log10(Total Bases)')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение Total Bases (log-scale)')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'sequencing_depth.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'sequencing_depth.png'}")
    
    def demographic_analysis(self):
        """Детальный анализ демографии"""
        print("\n" + "="*80)
        print("ДЕМОГРАФИЧЕСКИЙ АНАЛИЗ")
        print("="*80)
        
        # Пол
        sex_counts = self.df['sample_sex'].value_counts()
        sex_total = sex_counts.sum()
        print(f"\nПол (всего заполнено: {sex_total} из {len(self.df)}, {sex_total/len(self.df)*100:.1f}%):")
        for sex, count in sex_counts.items():
            pct = (count / len(self.df)) * 100
            print(f"  {sex:15s}: {count:6d} ({pct:5.1f}%)")
        
        # Возраст
        age_col = 'sample_age'
        age_data = pd.to_numeric(self.df[age_col], errors='coerce').dropna()
        print(f"\nВозраст (заполнено: {len(age_data)} из {len(self.df)}, {len(age_data)/len(self.df)*100:.1f}%):")
        if len(age_data) > 0:
            print(f"  Среднее:   {age_data.mean():.1f} лет")
            print(f"  Медиана:   {age_data.median():.1f} лет")
            print(f"  Мин:       {age_data.min():.0f} лет")
            print(f"  Макс:      {age_data.max():.0f} лет")
            print(f"  Ст. откл.: {age_data.std():.1f} лет")
        
        # Графики
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Sex distribution
        ax = axes[0, 0]
        if len(sex_counts) > 0:
            colors = sns.color_palette('Set2', n_colors=len(sex_counts))
            ax.pie(sex_counts.values, labels=sex_counts.index, autopct='%1.1f%%', 
                   startangle=90, colors=colors)
            ax.set_title('Распределение по полу (Illumina)')
        
        # Sex bar chart
        ax = axes[0, 1]
        if len(sex_counts) > 0:
            sex_counts.plot(kind='bar', ax=ax, color=colors)
            ax.set_xlabel('Пол')
            ax.set_ylabel('Количество образцов')
            ax.set_title('Распределение по полу')
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
            for i, v in enumerate(sex_counts.values):
                ax.text(i, v, f'{v}', ha='center', va='bottom')
        
        # Age distribution
        ax = axes[1, 0]
        if len(age_data) > 0:
            age_data.hist(bins=30, ax=ax, edgecolor='black')
            ax.set_xlabel('Возраст (лет)')
            ax.set_ylabel('Количество образцов')
            ax.set_title('Распределение по возрасту')
            ax.axvline(age_data.mean(), color='red', linestyle='--', 
                      label=f'Mean: {age_data.mean():.1f}')
            ax.axvline(age_data.median(), color='green', linestyle='--', 
                      label=f'Median: {age_data.median():.1f}')
            ax.legend()
        
        # Age boxplot by sex
        ax = axes[1, 1]
        if len(age_data) > 0 and len(sex_counts) > 0:
            age_sex_df = pd.DataFrame({
                'age': pd.to_numeric(self.df['sample_age'], errors='coerce'),
                'sex': self.df['sample_sex']
            }).dropna()
            
            if len(age_sex_df) > 0:
                age_sex_df.boxplot(column='age', by='sex', ax=ax)
                ax.set_xlabel('Пол')
                ax.set_ylabel('Возраст (лет)')
                ax.set_title('Распределение возраста по полу')
                plt.suptitle('')  # Remove automatic title
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'demographics.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'demographics.png'}")
    
    def disease_analysis(self):
        """Детальный анализ заболеваний"""
        print("\n" + "="*80)
        print("АНАЛИЗ ЗАБОЛЕВАНИЙ")
        print("="*80)
        
        # Ищем все колонки с информацией о заболеваниях
        disease_cols = [c for c in self.df.columns if 'disease' in c.lower()]
        print(f"\nНайдено колонок с информацией о заболеваниях: {len(disease_cols)}")
        
        # Основная колонка - sample_disease
        diseases = self.df['sample_disease'].value_counts()
        disease_total = diseases.sum()
        print(f"\nЗаболевания из sample_disease (заполнено: {disease_total} из {len(self.df)}, {disease_total/len(self.df)*100:.1f}%):")
        
        # Топ 20 заболеваний
        for disease, count in diseases.head(20).items():
            pct = (count / len(self.df)) * 100
            disease_str = str(disease)[:60] + '...' if len(str(disease)) > 60 else str(disease)
            print(f"  {disease_str:65s}: {count:6d} ({pct:5.1f}%)")
        
        if len(diseases) > 20:
            print(f"  ... и еще {len(diseases) - 20} заболеваний")
        
        # График топ заболеваний
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))
        
        # Топ 15 заболеваний - bar chart
        top_diseases = diseases.head(15)
        ax1.barh(range(len(top_diseases)), top_diseases.values)
        ax1.set_yticks(range(len(top_diseases)))
        labels = [str(d)[:50] + '...' if len(str(d)) > 50 else str(d) for d in top_diseases.index]
        ax1.set_yticklabels(labels, fontsize=9)
        ax1.set_xlabel('Количество образцов')
        ax1.set_title('Топ 15 заболеваний (Illumina)')
        ax1.invert_yaxis()
        
        for i, v in enumerate(top_diseases.values):
            ax1.text(v, i, f' {v}', va='center')
        
        # Pie chart для топ 10
        top10 = diseases.head(10)
        other = diseases[10:].sum()
        if other > 0:
            plot_data = pd.concat([top10, pd.Series({'Other': other})])
        else:
            plot_data = top10
        
        colors = sns.color_palette('tab20', n_colors=len(plot_data))
        ax2.pie(plot_data.values, labels=[str(l)[:30] for l in plot_data.index], 
                autopct='%1.1f%%', startangle=90, colors=colors)
        ax2.set_title('Доля основных заболеваний')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'diseases.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'diseases.png'}")
    
    def tissue_analysis(self):
        """Детальный анализ тканей"""
        print("\n" + "="*80)
        print("АНАЛИЗ ТКАНЕЙ")
        print("="*80)
        
        # sample_tissue и sample_body_site
        tissues = self.df['sample_tissue'].value_counts()
        tissue_total = tissues.sum()
        print(f"\nТкани из sample_tissue (заполнено: {tissue_total} из {len(self.df)}, {tissue_total/len(self.df)*100:.1f}%):")
        
        for tissue, count in tissues.head(20).items():
            pct = (count / len(self.df)) * 100
            tissue_str = str(tissue)[:60] + '...' if len(str(tissue)) > 60 else str(tissue)
            print(f"  {tissue_str:65s}: {count:6d} ({pct:5.1f}%)")
        
        if len(tissues) > 20:
            print(f"  ... и еще {len(tissues) - 20} типов тканей")
        
        # График
        fig, axes = plt.subplots(2, 1, figsize=(16, 12))
        
        # Топ 15 тканей
        top_tissues = tissues.head(15)
        ax = axes[0]
        ax.barh(range(len(top_tissues)), top_tissues.values)
        ax.set_yticks(range(len(top_tissues)))
        labels = [str(t)[:50] + '...' if len(str(t)) > 50 else str(t) for t in top_tissues.index]
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel('Количество образцов')
        ax.set_title('Топ 15 типов тканей (Illumina)')
        ax.invert_yaxis()
        
        for i, v in enumerate(top_tissues.values):
            ax.text(v, i, f' {v}', va='center')
        
        # Pie chart для топ 10
        ax = axes[1]
        top10 = tissues.head(10)
        other = tissues[10:].sum()
        if other > 0:
            plot_data = pd.concat([top10, pd.Series({'Other': other})])
        else:
            plot_data = top10
        
        colors = sns.color_palette('tab20', n_colors=len(plot_data))
        ax.pie(plot_data.values, labels=[str(l)[:30] for l in plot_data.index], 
               autopct='%1.1f%%', startangle=90, colors=colors)
        ax.set_title('Доля основных типов тканей')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'tissues.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'tissues.png'}")
    
    def completeness_analysis(self):
        """Анализ полноты данных"""
        print("\n" + "="*80)
        print("АНАЛИЗ ПОЛНОТЫ ДАННЫХ")
        print("="*80)
        
        # Подсчет заполненности
        self.df['completeness'] = self.df[self.sample_cols].notna().sum(axis=1)
        self.df['completeness_pct'] = (self.df['completeness'] / len(self.sample_cols)) * 100
        
        stats = self.df['completeness_pct'].describe()
        print(f"\nСтатистика по полноте данных:")
        print(f"  Среднее:     {stats['mean']:.2f}%")
        print(f"  Медиана:     {stats['50%']:.2f}%")
        print(f"  Мин:         {stats['min']:.2f}%")
        print(f"  Макс:        {stats['max']:.2f}%")
        print(f"  Ст. откл.:   {stats['std']:.2f}%")
        
        # Топ 10 самых полных образцов
        print(f"\nТоп 10 самых полных Illumina образцов:")
        top10 = self.df.nlargest(10, 'completeness')
        for i, (idx, row) in enumerate(top10.iterrows(), 1):
            print(f"  {i:2d}. {row['srr_id']:15s}: {int(row['completeness']):3d}/{len(self.sample_cols)} ({row['completeness_pct']:5.2f}%)")
        
        # График
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Histogram
        ax = axes[0, 0]
        self.df['completeness_pct'].hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel('Полнота данных (%)')
        ax.set_ylabel('Количество образцов')
        ax.set_title('Распределение полноты данных')
        ax.axvline(stats['mean'], color='red', linestyle='--', label=f"Mean: {stats['mean']:.1f}%")
        ax.axvline(stats['50%'], color='green', linestyle='--', label=f"Median: {stats['50%']:.1f}%")
        ax.legend()
        
        # Boxplot
        ax = axes[0, 1]
        self.df['completeness_pct'].plot(kind='box', ax=ax)
        ax.set_ylabel('Полнота данных (%)')
        ax.set_title('Boxplot полноты данных')
        
        # Полнота по инструментам
        ax = axes[1, 0]
        top_instruments = self.df['instrument_model'].value_counts().head(5).index
        completeness_by_instrument = self.df[self.df['instrument_model'].isin(top_instruments)].groupby('instrument_model')['completeness_pct'].mean().sort_values(ascending=False)
        completeness_by_instrument.plot(kind='barh', ax=ax)
        ax.set_xlabel('Средняя полнота данных (%)')
        ax.set_ylabel('Инструмент')
        ax.set_title('Полнота данных по топ-5 инструментам')
        
        # Полнота по типу библиотеки
        ax = axes[1, 1]
        completeness_by_strategy = self.df.groupby('library_strategy')['completeness_pct'].mean().sort_values(ascending=False).head(10)
        completeness_by_strategy.plot(kind='barh', ax=ax)
        ax.set_xlabel('Средняя полнота данных (%)')
        ax.set_ylabel('Стратегия библиотеки')
        ax.set_title('Полнота данных по стратегиям библиотек (топ-10)')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'completeness.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Сохранено: {self.output_dir / 'completeness.png'}")
    
    def temporal_analysis(self):
        """Временной анализ публикаций"""
        print("\n" + "="*80)
        print("ВРЕМЕННОЙ АНАЛИЗ")
        print("="*80)
        
        # Конвертируем даты
        self.df['run_published_date'] = pd.to_datetime(self.df['run_published'], errors='coerce')
        
        valid_dates = self.df['run_published_date'].dropna()
        print(f"\nДоступно дат публикации: {len(valid_dates)} из {len(self.df)} ({len(valid_dates)/len(self.df)*100:.1f}%)")
        
        if len(valid_dates) > 0:
            print(f"Период: {valid_dates.min().date()} - {valid_dates.max().date()}")
            
            # По годам
            self.df['year'] = self.df['run_published_date'].dt.year
            yearly = self.df['year'].value_counts().sort_index()
            
            print(f"\nОбразцов по годам:")
            for year, count in yearly.items():
                if pd.notna(year):
                    print(f"  {int(year)}: {count:6d}")
            
            # График
            fig, axes = plt.subplots(2, 1, figsize=(16, 10))
            
            # По годам
            ax = axes[0]
            yearly.plot(kind='bar', ax=ax)
            ax.set_xlabel('Год')
            ax.set_ylabel('Количество образцов')
            ax.set_title('Публикации Illumina образцов по годам')
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
            
            # По месяцам (последние 2 года)
            ax = axes[1]
            self.df['year_month'] = self.df['run_published_date'].dt.to_period('M')
            recent = self.df[self.df['run_published_date'] >= pd.Timestamp.now() - pd.DateOffset(years=2)]
            if len(recent) > 0:
                monthly = recent['year_month'].value_counts().sort_index()
                monthly.index = monthly.index.astype(str)
                monthly.plot(kind='line', ax=ax, marker='o')
                ax.set_xlabel('Месяц')
                ax.set_ylabel('Количество образцов')
                ax.set_title('Публикации за последние 2 года (помесячно)')
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'temporal.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Сохранено: {self.output_dir / 'temporal.png'}")
    
    def save_filtered_data(self):
        """Сохранение отфильтрованных Illumina данных"""
        output_file = self.output_dir / 'illumina_samples.csv'
        self.df.to_csv(output_file, index=False)
        print(f"\n✓ Illumina образцы сохранены: {output_file}")
        print(f"  Всего: {len(self.df)} образцов")
        print(f"  Размер: {output_file.stat().st_size / 1024 / 1024:.1f} MB")
    
    def generate_report(self):
        """Генерация текстового отчета"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        report_file = self.output_dir / f'illumina_eda_report_{timestamp}.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("ДЕТАЛЬНЫЙ EDA ОТЧЕТ - ILLUMINA ОБРАЗЦЫ\n")
            f.write("="*80 + "\n\n")
            f.write(f"Дата создания: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Исходный файл: {len(self.df_full)} образцов (всего)\n")
            f.write(f"Illumina образцов: {len(self.df)}\n\n")
            
            # Платформы
            f.write("="*80 + "\n")
            f.write("ПЛАТФОРМЫ ILLUMINA\n")
            f.write("="*80 + "\n")
            instruments = self.df['instrument_model'].value_counts()
            for instrument, count in instruments.items():
                pct = (count / len(self.df)) * 100
                f.write(f"{instrument:40s}: {count:6d} ({pct:5.1f}%)\n")
            
            # Стратегии библиотек
            f.write("\n" + "="*80 + "\n")
            f.write("СТРАТЕГИИ БИБЛИОТЕК\n")
            f.write("="*80 + "\n")
            strategies = self.df['library_strategy'].value_counts()
            for strategy, count in strategies.items():
                pct = (count / len(self.df)) * 100
                f.write(f"{strategy:30s}: {count:6d} ({pct:5.1f}%)\n")
            
            # Демография
            f.write("\n" + "="*80 + "\n")
            f.write("ДЕМОГРАФИЯ\n")
            f.write("="*80 + "\n")
            sex_counts = self.df['sample_sex'].value_counts()
            f.write(f"Пол (заполнено: {sex_counts.sum()}/{len(self.df)}):\n")
            for sex, count in sex_counts.items():
                pct = (count / len(self.df)) * 100
                f.write(f"  {sex:15s}: {count:6d} ({pct:5.1f}%)\n")
            
            age_data = pd.to_numeric(self.df['sample_age'], errors='coerce').dropna()
            f.write(f"\nВозраст (заполнено: {len(age_data)}/{len(self.df)}):\n")
            if len(age_data) > 0:
                f.write(f"  Среднее:   {age_data.mean():.1f} лет\n")
                f.write(f"  Медиана:   {age_data.median():.1f} лет\n")
                f.write(f"  Диапазон:  {age_data.min():.0f} - {age_data.max():.0f} лет\n")
            
            # Заболевания
            f.write("\n" + "="*80 + "\n")
            f.write("ЗАБОЛЕВАНИЯ (ТОП 20)\n")
            f.write("="*80 + "\n")
            diseases = self.df['sample_disease'].value_counts()
            for disease, count in diseases.head(20).items():
                pct = (count / len(self.df)) * 100
                disease_str = str(disease)[:60]
                f.write(f"{disease_str:65s}: {count:6d} ({pct:5.1f}%)\n")
            
            # Ткани
            f.write("\n" + "="*80 + "\n")
            f.write("ТКАНИ (ТОП 20)\n")
            f.write("="*80 + "\n")
            tissues = self.df['sample_tissue'].value_counts()
            for tissue, count in tissues.head(20).items():
                pct = (count / len(self.df)) * 100
                tissue_str = str(tissue)[:60]
                f.write(f"{tissue_str:65s}: {count:6d} ({pct:5.1f}%)\n")
            
            # Полнота данных
            f.write("\n" + "="*80 + "\n")
            f.write("ПОЛНОТА ДАННЫХ\n")
            f.write("="*80 + "\n")
            stats = self.df['completeness_pct'].describe()
            f.write(f"Среднее:     {stats['mean']:.2f}%\n")
            f.write(f"Медиана:     {stats['50%']:.2f}%\n")
            f.write(f"Мин:         {stats['min']:.2f}%\n")
            f.write(f"Макс:        {stats['max']:.2f}%\n")
            
            f.write("\nТоп 10 самых полных образцов:\n")
            top10 = self.df.nlargest(10, 'completeness')
            for i, (idx, row) in enumerate(top10.iterrows(), 1):
                f.write(f"  {i:2d}. {row['srr_id']:15s}: {int(row['completeness']):3d}/{len(self.sample_cols)} ({row['completeness_pct']:5.2f}%)\n")
        
        print(f"\n✓ Отчет сохранен: {report_file}")
    
    def run_full_analysis(self):
        """Запуск полного анализа"""
        print("\n" + "="*80)
        print("ЗАПУСК ПОЛНОГО EDA ДЛЯ ILLUMINA ОБРАЗЦОВ")
        print("="*80 + "\n")
        
        self.platform_overview()
        self.library_strategy_analysis()
        self.library_layout_analysis()
        self.sequencing_depth_analysis()
        self.demographic_analysis()
        self.disease_analysis()
        self.tissue_analysis()
        self.completeness_analysis()
        self.temporal_analysis()
        self.save_filtered_data()
        self.generate_report()
        
        print("\n" + "="*80)
        print("АНАЛИЗ ЗАВЕРШЕН!")
        print("="*80)
        print(f"\nРезультаты сохранены в: {self.output_dir}")
        print(f"\nСоздано файлов:")
        for file in sorted(self.output_dir.glob('*')):
            size = file.stat().st_size
            if size > 1024*1024:
                size_str = f"{size/1024/1024:.1f} MB"
            elif size > 1024:
                size_str = f"{size/1024:.1f} KB"
            else:
                size_str = f"{size} B"
            print(f"  - {file.name:50s} ({size_str})")


def main():
    """Главная функция"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Детальный EDA для Illumina образцов')
    parser.add_argument('--metadata', type=str, 
                       default='data/metadata/sra_metadata_complete_20251117_085704.csv',
                       help='Путь к файлу с метаданными')
    
    args = parser.parse_args()
    
    # Создаем и запускаем анализ
    eda = IlluminaEDA(args.metadata)
    eda.run_full_analysis()


if __name__ == '__main__':
    main()
