#!/bin/bash
# Quick view of current metadata download progress

PARTIAL_FILE="/home/nicolaedrabcinski/cardiogen/data/metadata/sra_metadata_PARTIAL_IN_PROGRESS.json"

if [ -f "$PARTIAL_FILE" ]; then
    COUNT=$(python3 -c "import json; data=json.load(open('$PARTIAL_FILE')); print(len(data))")
    SIZE=$(du -h "$PARTIAL_FILE" | cut -f1)
    
    echo "========================================"
    echo "ТЕКУЩИЙ ПРОГРЕСС ЗАГРУЗКИ МЕТАДАННЫХ"
    echo "========================================"
    echo "Файл: sra_metadata_PARTIAL_IN_PROGRESS.json"
    echo "Сохранено образцов: $COUNT / 4596"
    echo "Размер файла: $SIZE"
    echo ""
    echo "Для просмотра данных:"
    echo "  python scripts/view_partial_metadata.py"
    echo ""
    echo "Для EDA (создаст графики):"
    echo "  python scripts/metadata_eda.py $PARTIAL_FILE"
    echo "========================================"
else
    echo "Файл $PARTIAL_FILE не найден"
    echo "Загрузка еще не началась или уже завершена"
fi
