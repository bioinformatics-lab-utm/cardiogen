#!/bin/bash
# Monitor metadata download progress

LOG_FILE="/home/nicolaedrabcinski/cardiogen/metadata_download.log"
PID_PATTERN="download_all_metadata.py"

echo "=========================================="
echo "METADATA DOWNLOAD MONITORING"
echo "=========================================="
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check if process is running
if pgrep -f "$PID_PATTERN" > /dev/null; then
    PID=$(pgrep -f "$PID_PATTERN")
    echo "Status: RUNNING (PID: $PID)"
    
    # Get process info
    ps -p $PID -o pid,ppid,%cpu,%mem,etime,cmd --no-headers
    echo ""
else
    echo "Status: NOT RUNNING"
    echo ""
fi

# Check log file
if [ -f "$LOG_FILE" ]; then
    echo "Log file: $LOG_FILE"
    echo "Log size: $(du -h $LOG_FILE | cut -f1)"
    echo ""
    
    # Extract progress
    TOTAL=$(grep -oP "Total unique SRR IDs: \K\d+" "$LOG_FILE" | tail -1)
    COMPLETED=$(grep -oP "Progress: \K\d+(?=/)" "$LOG_FILE" | tail -1)
    
    if [ ! -z "$TOTAL" ] && [ ! -z "$COMPLETED" ]; then
        PERCENT=$(awk "BEGIN {printf \"%.2f\", ($COMPLETED/$TOTAL)*100}")
        echo "Progress: $COMPLETED / $TOTAL SRR IDs ($PERCENT%)"
        
        # Estimate time remaining
        if [ -f "$LOG_FILE" ]; then
            START_TIME=$(stat -c %Y "$LOG_FILE")
            CURRENT_TIME=$(date +%s)
            ELAPSED=$((CURRENT_TIME - START_TIME))
            
            if [ $COMPLETED -gt 0 ]; then
                TIME_PER_ITEM=$(awk "BEGIN {printf \"%.2f\", $ELAPSED/$COMPLETED}")
                REMAINING=$((TOTAL - COMPLETED))
                ETA_SECONDS=$(awk "BEGIN {printf \"%.0f\", $TIME_PER_ITEM * $REMAINING}")
                ETA_HOURS=$(awk "BEGIN {printf \"%.1f\", $ETA_SECONDS/3600}")
                
                echo "Elapsed time: $(($ELAPSED / 3600))h $(($ELAPSED % 3600 / 60))m"
                echo "Estimated time remaining: ${ETA_HOURS}h"
                echo "Average: ${TIME_PER_ITEM}s per SRR ID"
            fi
        fi
    else
        echo "Progress: Initializing..."
    fi
    echo ""
    
    # Show recent errors
    ERROR_COUNT=$(grep -c "429 Client Error" "$LOG_FILE" 2>/dev/null || echo "0")
    echo "Rate limit errors (429): $ERROR_COUNT"
    echo ""
    
    # Show last few lines
    echo "Last 10 lines of log:"
    echo "---"
    tail -10 "$LOG_FILE"
else
    echo "Log file not found: $LOG_FILE"
fi

echo ""
echo "=========================================="

# Check metadata output directory
METADATA_DIR="/home/nicolaedrabcinski/cardiogen/data/metadata"
if [ -d "$METADATA_DIR" ]; then
    echo "Metadata directory: $METADATA_DIR"
    echo "Files:"
    ls -lh "$METADATA_DIR" 2>/dev/null | tail -n +2
fi

echo "=========================================="
