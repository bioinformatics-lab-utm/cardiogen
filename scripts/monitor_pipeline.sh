#!/bin/bash
# Мониторинг запущенного пайплайна

echo "🔬 Cardiogen Pipeline Monitoring"
echo "================================"

# Проверка процесса
if [ -f pipeline.pid ]; then
    PID=$(cat pipeline.pid)
    if ps -p $PID > /dev/null 2>&1; then
        echo "✅ Pipeline is RUNNING (PID: $PID)"
        
        # Время работы
        ELAPSED=$(ps -p $PID -o etime= | tr -d ' ')
        echo "⏱️  Running time: $ELAPSED"
        
        # Использование CPU и памяти
        CPU=$(ps -p $PID -o %cpu= | tr -d ' ')
        MEM=$(ps -p $PID -o %mem= | tr -d ' ')
        echo "💻 CPU: ${CPU}% | Memory: ${MEM}%"
    else
        echo "❌ Pipeline STOPPED (PID $PID not found)"
        echo "   Check pipeline_full.log for details"
    fi
else
    echo "⚠️  No pipeline.pid found"
    echo "   Is the pipeline running?"
fi

echo ""
echo "📊 Progress Statistics"
echo "================================"

# Подсчёт обработанных батчей
BATCHES=$(ls -1 data/results/batch_*_results.json 2>/dev/null | wc -l)
echo "Batches completed: $BATCHES"

# Подсчёт QC отчётов
QC_REPORTS=$(find data/qc -name "*_fastqc.html" 2>/dev/null | wc -l)
echo "QC reports generated: $QC_REPORTS"

# Использование диска
echo ""
echo "💾 Disk Usage"
echo "================================"
echo "QC directory: $(du -sh data/qc 2>/dev/null | cut -f1)"
echo "Results directory: $(du -sh data/results 2>/dev/null | cut -f1)"
echo "Total project: $(du -sh . 2>/dev/null | cut -f1)"

# Последние строки лога
echo ""
echo "📋 Last 10 log lines"
echo "================================"
tail -10 pipeline_full.log 2>/dev/null || echo "No log file found"

echo ""
echo "Commands:"
echo "  tail -f pipeline_full.log    # Watch live logs"
echo "  ./scripts/monitor_pipeline.sh  # Re-run this monitor"
echo "  cat data/results/summary_report.txt  # View summary"
