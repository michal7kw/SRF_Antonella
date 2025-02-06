#!/bin/bash
#SBATCH --job-name=5_differential_binding_deseq2
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding_deseq2.err"
#SBATCH --output="logs/5_differential_binding_deseq2.out"

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to validate and clean peak file format
validate_peak_file() {
    local peak_file=$1
    local expected_columns=9
    local temp_file="${peak_file}.tmp"
    
    # Check if file exists
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Peak file not found: $peak_file"
        return 1
    fi
    
    # Clean and validate the file
    awk -v cols=$expected_columns '
    BEGIN { valid_lines = 0; total_lines = 0 }
    {
        total_lines++
        # Skip empty lines or lines with wrong number of columns
        if (NF == cols) {
            print $0
            valid_lines++
        }
    }
    END {
        if (valid_lines == 0) {
            printf "ERROR: No valid lines found in file\n" > "/dev/stderr"
            exit 1
        }
        if (valid_lines < total_lines) {
            printf "WARNING: Filtered %d invalid lines from %d total lines\n", \
                   total_lines - valid_lines, total_lines > "/dev/stderr"
        }
    }' "$peak_file" > "$temp_file" || {
        log_message "ERROR: Peak file validation failed for $peak_file"
        rm -f "$temp_file"
        return 1
    }
    
    # Replace original with cleaned file
    mv "$temp_file" "$peak_file"
    return 0
}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

PEAKS_SUFFIX="peaks_final"
# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs
mkdir -p analysis/diffbind_broad_deseq2
mkdir -p analysis/annotation_broad_deseq2

# Define sample list with only selected samples that have similar peak counts
samples=(GFP_1 GFP_3 YAF_2 YAF_3)  # Samples with similar peak counts (27k-37k peaks)
type="broad"

log_message "Using selected samples with similar peak counts:"

# Validate all input files first
log_message "Validating input files for all samples..."
for sample in "${samples[@]}"; do
    peak_file="analysis/peaks/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
    bam_file="analysis/aligned/${sample}.dedup.bam"
    
    if [[ ! -f "$bam_file" ]]; then
        log_message "ERROR: BAM file not found: $bam_file"
        exit 1
    fi
    
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Peak file not found: $peak_file"
        exit 1
    fi
    
    # Create a backup of the peak file
    peak_file_backup="${peak_file}.backup"
    if ! cp "$peak_file" "$peak_file_backup"; then
        log_message "ERROR: Failed to create backup of peak file"
        exit 1
    fi
    log_message "Created backup of peak file for ${sample}"
done

# Validate and clean all peak files
for sample in "${samples[@]}"; do
    peak_file="analysis/peaks/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
    peak_file_backup="${peak_file}.backup"
    
    log_message "Validating peak file format for ${sample}..."
    if ! validate_peak_file "$peak_file"; then
        log_message "ERROR: Peak file validation failed. Please check the format of $peak_file"
        # Restore from backup
        mv "$peak_file_backup" "$peak_file"
        exit 1
    fi
    
    # If successful, remove backup
    rm -f "$peak_file_backup"
    log_message "Peak validation completed for ${sample}"
done

# Run R script for differential binding analysis
log_message "Running differential binding analysis..."
if ! Rscript scripts/5_differential_binding_deseq2.R; then
    log_message "ERROR: Differential binding analysis failed"
    exit 1
fi

log_message "Differential binding analysis completed successfully" 