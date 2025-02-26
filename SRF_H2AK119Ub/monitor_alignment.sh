#!/bin/bash
# Script to monitor alignment jobs and detect potential hanging

# Set the maximum time (in seconds) an alignment job should take
MAX_RUNTIME=7200  # 2 hours

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo "=== Alignment Job Monitor ==="
echo "Checking for alignment jobs every 5 minutes"
echo "Press Ctrl+C to exit"
echo

while true; do
    # Get current timestamp
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
    
    # Get all running alignment jobs
    ALIGN_JOBS=$(squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R" | grep -i "align")
    
    if [ -z "$ALIGN_JOBS" ]; then
        echo -e "${TIMESTAMP} - ${GREEN}No alignment jobs currently running${NC}"
    else
        echo -e "${TIMESTAMP} - ${YELLOW}Found running alignment jobs:${NC}"
        echo "$ALIGN_JOBS"
        
        # Check log files for potential issues
        echo -e "\nChecking alignment log files for potential issues:"
        
        # Find all alignment log files modified in the last 24 hours
        LOG_FILES=$(find logs/alignment -name "*.align.log" -mtime -1 2>/dev/null)
        
        for LOG in $LOG_FILES; do
            SAMPLE=$(basename "$LOG" .align.log)
            
            # Check if the log file has been updated in the last 15 minutes
            if [ $(find "$LOG" -mmin -15 -print 2>/dev/null | wc -l) -eq 0 ]; then
                # Check if the log file contains a completion message
                if grep -q "Alignment rule completed successfully" "$LOG"; then
                    echo -e "${GREEN}$SAMPLE: Completed successfully${NC}"
                else
                    # Check the last modification time
                    LAST_MOD=$(stat -c %Y "$LOG")
                    CURRENT=$(date +%s)
                    ELAPSED=$((CURRENT - LAST_MOD))
                    
                    if [ $ELAPSED -gt $MAX_RUNTIME ]; then
                        echo -e "${RED}$SAMPLE: WARNING - Log file hasn't been updated in $(($ELAPSED/60)) minutes, job might be hanging${NC}"
                        echo -e "${YELLOW}Last 5 lines of $LOG:${NC}"
                        tail -n 5 "$LOG"
                    else
                        echo -e "${YELLOW}$SAMPLE: In progress - Last updated $(($ELAPSED/60)) minutes ago${NC}"
                    fi
                fi
            else
                echo -e "${GREEN}$SAMPLE: Active - Log updated recently${NC}"
            fi
        done
    fi
    
    echo -e "\nSleeping for 5 minutes before next check..."
    sleep 300
done
