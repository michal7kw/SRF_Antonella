#!/bin/bash

# Documentation:
# This script performs differential binding analysis for CUT&Tag data using DiffBind
# It processes multiple samples with similar peak counts to identify differentially bound regions
# The script performs the following main steps:
# 1. Validates input peak and BAM files
# 2. Cleans and standardizes peak file formats
# 3. Runs differential binding analysis using an R script
# 4. Generates output files with differential binding results
# This version is adapted for local execution.

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
# Validate a broadPeak file format and clean it to ensure consistency for DiffBind analysis
# The function performs the following checks:
# 1. Checks if the file exists
# 2. Checks for empty lines
# 3. Checks the number of columns (should be 9)
# 4. Validates the chromosome field (accepting both "chr1" and "1" formats)
# 5. Validates the numeric fields (start, end, score, signalValue, pValue, qValue)
# 6. Validates the strand field
# If any errors are found, the function prints an error message and returns 1
# Otherwise, it returns 0 and replaces the original file with a cleaned version

validate_peak_file() {
    local peak_file=$1
    local expected_columns=9  # Standard number of columns in broadPeak format
    local temp_file="${peak_file}.tmp"
    local error_file="${peak_file}.errors"

    # Check if file exists
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Peak file not found: $peak_file"
        return 1
    fi

    # Clean and validate the file using awk with detailed error reporting
    awk -v cols=$expected_columns -v error_file="$error_file" '
    BEGIN {
        valid_lines = 0;
        total_lines = 0;
        print "Error report for " FILENAME > error_file
    }
    {
        total_lines++
        # Skip empty lines
        if (NF == 0) {
            printf "Line %d: Empty line\n", NR > error_file
            next
        }

        # Check number of columns
        if (NF != cols) {
            printf "Line %d: Found %d columns, expected %d\n", NR, NF, cols > error_file
            next
        }

        # Validate chromosome field (accepting both "chr1" and "1" formats)
        chr_value = $1
        # Remove "chr" prefix if present for validation
        sub(/^chr/, "", chr_value)
        if (chr_value !~ /^([1-9][0-9]?|X|Y|M|MT)$/) {
            printf "Line %d: Invalid chromosome format: %s\n", NR, $1 > error_file
            next
        }

        # Validate numeric fields (start, end, score, signalValue, pValue, qValue)
        if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ ||
            $7 !~ /^[0-9.]+$/ || $8 !~ /^[0-9.]+$/ || $9 !~ /^[0-9.]+$/) {
            printf "Line %d: Invalid numeric field(s)\n", NR > error_file
            next
        }

        # Validate strand field
        if ($6 != "." && $6 != "+" && $6 != "-") {
            printf "Line %d: Invalid strand: %s\n", NR, $6 > error_file
            next
        }

        # If all validations pass, print the line
        print $0
        valid_lines++
    }
    END {
        printf "\nSummary:\n" > error_file
        printf "Total lines processed: %d\n", total_lines > error_file
        printf "Valid lines: %d\n", valid_lines > error_file
        printf "Invalid lines: %d\n", total_lines - valid_lines > error_file

        if (valid_lines == 0) {
            printf "ERROR: No valid lines found in file\n" > "/dev/stderr"
            exit 1
        }
    }' "$peak_file" > "$temp_file" || {
        log_message "ERROR: Peak file validation failed for $peak_file"
        log_message "Check ${error_file} for detailed error report"
        rm -f "$temp_file"
        return 1
    }

    # If validation successful, replace original with cleaned file
    mv "$temp_file" "$peak_file"

    # Report validation results
    log_message "Peak file validation completed for $peak_file"
    log_message "See ${error_file} for validation report"
    return 0
}

# Activate conda environment with required tools
# Ensure the 'snakemake' conda environment is active before running this script
# Example: conda activate snakemake
# source /opt/common/tools/ric.cosr/miniconda3/bin/activate # Removed cluster-specific path
# conda activate snakemake # Assuming environment is already active

# Define constants
PEAKS_SUFFIX="peaks_final"  # Suffix for final peak files
type="broad"  # Analysis type (broad peaks)

# Define working directory (assuming script is run from 1_iterative_processing)
WORKDIR="."
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define directories
OUTPUT_DIR="analysis/7_differential_binding_v2"
PEAKS_DIR="analysis/5_peak_calling_v2"
ALIGN_DIR="analysis/3_alignment"
LOG_DIR="logs"

# Create necessary output directories
log_message "Creating output directories..."
mkdir -p ${LOG_DIR}  # For log files
mkdir -p ${OUTPUT_DIR}  # For DiffBind results

# Define sample information
SAMPLE_IDS="GFP_1,GFP_2,GFP_3,YAF_1,YAF_2,YAF_3"
SAMPLE_CONDITIONS="GFP,GFP,GFP,YAF,YAF,YAF"
SAMPLE_REPLICATES="1,2,3,1,2,3"

# Define sample list with only selected samples that have similar peak counts
# Using samples with comparable peak counts (27k-37k peaks) for more reliable comparison
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

log_message "Using selected samples with similar peak counts: ${samples[*]}"

# Validate all input files first
log_message "Validating input files for all samples..."
for sample in "${samples[@]}"; do
    peak_file="${PEAKS_DIR}/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
    bam_file="${ALIGN_DIR}/${sample}.dedup.bam"

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
        log_message "ERROR: Failed to create backup of peak file for $sample"
        exit 1
    fi
    log_message "Created backup of peak file for ${sample}"
done

# Validate and clean all peak files
for sample in "${samples[@]}"; do
    peak_file="${PEAKS_DIR}/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
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
# Redirect stdout and stderr to log files
if ! Rscript scripts/7_differential_binding_v2.R \
    "${OUTPUT_DIR}" \
    "${PEAKS_DIR}" \
    "${PEAKS_SUFFIX}" \
    "${ALIGN_DIR}" \
    "${SAMPLE_IDS}" \
    "${SAMPLE_CONDITIONS}" \
    "${SAMPLE_REPLICATES}" > "${LOG_DIR}/7_differential_binding.out" 2> "${LOG_DIR}/7_differential_binding.err"; then
    log_message "ERROR: Differential binding analysis failed. Check logs/${LOG_DIR}/7_differential_binding.err"
    exit 1
fi

log_message "Differential binding analysis completed successfully"
log_message "Output files are in ${OUTPUT_DIR}"
log_message "Logs are in ${LOG_DIR}"