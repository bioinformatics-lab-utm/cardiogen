#!/bin/bash
# Cleanup script for cardiogen project

echo "🧹 Cardiogen Cleanup Script"
echo "=========================="

# Function to get directory size
get_size() {
    du -sh "$1" 2>/dev/null | cut -f1
}

echo ""
echo "Current disk usage:"
df -h . | grep -v "Filesystem"

echo ""
echo "Directory sizes:"
echo "  archive/: $(get_size archive)"
echo "  data/: $(get_size data)"
echo "  data/qc/: $(get_size data/qc)"
echo "  data/results/: $(get_size data/results)"

echo ""
read -p "What would you like to clean? (1=logs, 2=temp files, 3=archive, 4=all, 0=cancel): " choice

case $choice in
    1)
        echo "Cleaning log files..."
        rm -f *.log
        echo "✅ Logs cleaned"
        ;;
    2)
        echo "Cleaning temporary files..."
        rm -rf data/temp_batch_*
        rm -rf fasterq.tmp.*
        rm -rf __pycache__
        echo "✅ Temporary files cleaned"
        ;;
    3)
        echo "⚠️  This will delete archived data!"
        read -p "Are you sure? (yes/no): " confirm
        if [ "$confirm" = "yes" ]; then
            rm -rf archive/
            echo "✅ Archive cleaned"
        else
            echo "❌ Cancelled"
        fi
        ;;
    4)
        echo "⚠️  This will delete logs, temp files, and archive!"
        read -p "Are you sure? (yes/no): " confirm
        if [ "$confirm" = "yes" ]; then
            rm -f *.log
            rm -rf data/temp_batch_*
            rm -rf fasterq.tmp.*
            rm -rf __pycache__
            rm -rf archive/
            echo "✅ Full cleanup completed"
        else
            echo "❌ Cancelled"
        fi
        ;;
    0)
        echo "❌ Cleanup cancelled"
        ;;
    *)
        echo "❌ Invalid choice"
        ;;
esac

echo ""
echo "Final disk usage:"
df -h . | grep -v "Filesystem"
