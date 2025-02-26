#!/bin/bash
#SBATCH --job-name=7_differential_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_differential_binding.err"
#SBATCH --output="logs/7_differential_binding.out"

# Documentation:
# This script performs differential binding analysis for CUT&Tag data using DiffBind
# It processes multiple samples with similar peak counts to identify differentially bound regions
# The script performs the following main steps:
# 1. Validates input peak and BAM files
# 2. Cleans and standardizes peak file formats
# 3. Runs differential binding analysis using an R script
# 4. Generates output files with differential binding results

# Input files:
# - analysis/5_peak_calling/{sample}_broad_peaks_final.broadPeak: Final peak calls for each sample
# - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample
# - scripts/7_differential_binding.R: R script containing DiffBind analysis code

# Output files:
# - analysis/7_differential_binding/: Directory containing DiffBind results
# - analysis/7_differential_binding/: Directory containing annotated differential binding results
# - logs/7_differential_binding.out: Standard output log
# - logs/7_differential_binding.err: Error log

# Error handling
set -e  # Exit immediately if a command exits with non-zero status
set -u  # Treat unset variables as an error
set -o pipefail  # Fail pipeline if any command fails

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to validate and clean peak file format
validate_peak_file() {
    local peak_file=$1
    local expected_columns=9  # Standard number of columns in broadPeak format
    local temp_file="${peak_file}.tmp"
    
    # Check if file exists
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Peak file not found: $peak_file"
        return 1
    fi
    
    # Clean and validate the file using awk
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

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define constants
PEAKS_SUFFIX="peaks_final"  # Suffix for final peak files
type="broad"  # Analysis type (broad peaks)

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/7_differential_binding"

# Create necessary output directories
log_message "Creating output directories..."
mkdir -p logs  # For log files
mkdir -p ${OUTPUT_DIR}  # For DiffBind results

# Define sample list with only selected samples that have similar peak counts
# Using samples with comparable peak counts (27k-37k peaks) for more reliable comparison
samples=(GFP_1 GFP_3 YAF_2 YAF_3)

log_message "Using selected samples with similar peak counts: ${samples[*]}"

# Validate all input files first
log_message "Validating input files for all samples..."
for sample in "${samples[@]}"; do
    peak_file="analysis/5_peak_calling/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
    bam_file="analysis/3_alignment/${sample}.dedup.bam"
    
    # Check if BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        log_message "ERROR: BAM file not found: $bam_file"
        exit 1
    fi
    
    # Check if peak file exists
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Peak file not found: $peak_file"
        exit 1
    fi
    
    # Create a backup of the peak file before validation
    peak_file_backup="${peak_file}.backup"
    if ! cp "$peak_file" "$peak_file_backup"; then
        log_message "ERROR: Failed to create backup of peak file"
        exit 1
    fi
    log_message "Created backup of peak file for ${sample}"
done

# Validate and clean all peak files
for sample in "${samples[@]}"; do
    peak_file="analysis/5_peak_calling/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
    peak_file_backup="${peak_file}.backup"
    
    log_message "Validating peak file format for ${sample}..."
    if ! validate_peak_file "$peak_file"; then
        log_message "ERROR: Peak file validation failed. Please check the format of $peak_file"
        # Restore from backup if validation fails
        mv "$peak_file_backup" "$peak_file"
        exit 1
    fi
    
    # If validation is successful, remove the backup
    rm -f "$peak_file_backup"
    log_message "Peak validation completed for ${sample}"
done

# Run R script for differential binding analysis
log_message "Running differential binding analysis..."
if ! Rscript scripts/7_differential_binding.R ${OUTPUT_DIR}; then
    log_message "ERROR: Differential binding analysis failed"
    exit 1
fi

log_message "Differential binding analysis completed successfully"