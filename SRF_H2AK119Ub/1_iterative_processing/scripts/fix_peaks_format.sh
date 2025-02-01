#!/bin/bash
#SBATCH --job-name=fix_peaks_format
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/fix_peaks_format.err"
#SBATCH --output="logs/fix_peaks_format.out"

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
PEAKS_DIR="$BASE_DIR/analysis/peaks"

log_message "Starting peak file format fix for all peak files"

# Process all peak files (both broad and narrow)
for FILE_TO_FIX in $PEAKS_DIR/*Peak; do
    log_message "Processing file: $(basename $FILE_TO_FIX)"
    
    # Create a backup first
    log_message "Creating backup of original file..."
    if cp "$FILE_TO_FIX" "$FILE_TO_FIX.backup"; then
        log_message "Backup created successfully"
    else
        log_message "ERROR: Failed to create backup for $FILE_TO_FIX"
        continue
    fi

    # Use awk to convert scientific notation to regular integers
    log_message "Converting scientific notation to integers..."
    if awk 'BEGIN{OFS="\t"} {
        for(i=1;i<=NF;i++) {
            if($i ~ /^[0-9]+\.[0-9]+e[+-][0-9]+$/) {
                $i = sprintf("%.0f", $i)
            }
        }
        print
    }' "$FILE_TO_FIX" > "$FILE_TO_FIX.tmp" && mv "$FILE_TO_FIX.tmp" "$FILE_TO_FIX"; then
        log_message "Successfully converted file format for $FILE_TO_FIX"
    else
        log_message "ERROR: Failed to convert file format for $FILE_TO_FIX"
        continue
    fi
done

log_message "Script completed successfully"
