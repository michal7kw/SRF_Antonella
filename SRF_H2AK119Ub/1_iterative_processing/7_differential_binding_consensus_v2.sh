#!/bin/bash

# Documentation:
# This script performs differential binding analysis for CUT&Tag data using DiffBind,
# focusing on consensus peaks (peaks present in at least two replicates per condition).
# It processes multiple samples to identify differentially bound regions based on these consensus peaks.
# The script performs the following main steps:
# 1. Validates input consensus peak and BAM files.
# 2. Runs differential binding analysis using an R script (7_differential_binding_consensus_v2.R).
# 3. Generates output files with differential binding results.
# This version is adapted for local execution and uses consensus peaks.

# Input files:
# - analysis/6_consensus_peaks/{condition}_consensus_peaks.bed: Consensus peak calls for each condition.
# - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample.
# - scripts/7_differential_binding_consensus_v2.R: R script containing DiffBind analysis code for consensus peaks.

# Output files:
# - analysis/7_differential_binding_consensus_v2/: Directory containing DiffBind results.
# - logs/7_differential_binding_consensus_v2.out: Standard output log.
# - logs/7_differential_binding_consensus_v2.err: Error log.

# Error handling
set -e  # Exit immediately if a command exits with non-zero status
set -u  # Treat unset variables as an error
set -o pipefail  # Fail pipeline if any command fails

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to validate a BED file (checks existence and non-emptiness)
validate_consensus_peak_file() {
    local peak_file=$1
    log_message "Validating consensus peak file: $peak_file"
    if [[ ! -f "$peak_file" ]]; then
        log_message "ERROR: Consensus peak file not found: $peak_file"
        return 1
    fi
    if [[ ! -s "$peak_file" ]]; then
        log_message "ERROR: Consensus peak file is empty: $peak_file"
        return 1
    fi
    # Basic check for at least 3 columns (chr, start, end)
    if awk 'NF < 3 {exit 1}' "$peak_file"; then
        : # File has at least 3 columns on all lines
    else
        log_message "ERROR: Consensus peak file $peak_file has lines with fewer than 3 columns."
        return 1
    fi
    log_message "Consensus peak file $peak_file appears valid."
    return 0
}

# Activate conda environment with required tools (if necessary)
# Ensure the environment containing R, DiffBind, and its dependencies is active.
# Example: conda activate snakemake_or_diffbind_env

# Define working directory (assuming script is run from 1_iterative_processing)
WORKDIR="."
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define directories
OUTPUT_DIR="analysis/7_differential_binding_consensus_v2"
CONSENSUS_PEAKS_DIR="analysis/6_consensus_peaks" # Directory where consensus peaks are stored
ALIGN_DIR="analysis/3_alignment"
LOG_DIR="logs"
SCRIPT_DIR="scripts"

# Create necessary output directories
log_message "Creating output directories..."
mkdir -p ${LOG_DIR}
mkdir -p ${OUTPUT_DIR}

# Define sample information for selected samples
# Using samples with comparable peak counts (example from v2 script)
# This list determines which BAM files are used.
declare -a samples_to_process=("GFP_1" "GFP_3" "YAF_2" "YAF_3")
# samples_to_process=("GFP_1" "GFP_2" "GFP_3" "YAF_1" "YAF_2" "YAF_3") # Uncomment to use all samples

log_message "Using selected samples for BAM files: ${samples_to_process[*]}"

# Define full sample metadata (used to derive info for the selected samples)
declare -A ALL_SAMPLES_INFO
ALL_SAMPLES_INFO[GFP_1]="GFP,1"
ALL_SAMPLES_INFO[GFP_2]="GFP,2"
ALL_SAMPLES_INFO[GFP_3]="GFP,3"
ALL_SAMPLES_INFO[YAF_1]="YAF,1"
ALL_SAMPLES_INFO[YAF_2]="YAF,2"
ALL_SAMPLES_INFO[YAF_3]="YAF,3"

# Dynamically build R_SAMPLE_IDS, R_SAMPLE_CONDITIONS, R_SAMPLE_REPLICATES for the R script
R_SAMPLE_IDS_ARR=()
R_SAMPLE_CONDITIONS_ARR=()
R_SAMPLE_REPLICATES_ARR=()

for sample_id in "${samples_to_process[@]}"; do
    if [[ -z "${ALL_SAMPLES_INFO[$sample_id]+_}" ]]; then
        log_message "ERROR: Sample ID $sample_id not found in ALL_SAMPLES_INFO. Please define it."
        exit 1
    fi
    R_SAMPLE_IDS_ARR+=("$sample_id")
    condition=$(echo "${ALL_SAMPLES_INFO[$sample_id]}" | cut -d',' -f1)
    replicate=$(echo "${ALL_SAMPLES_INFO[$sample_id]}" | cut -d',' -f2)
    R_SAMPLE_CONDITIONS_ARR+=("$condition")
    R_SAMPLE_REPLICATES_ARR+=("$replicate")
done

R_PASSED_SAMPLE_IDS=$(IFS=,; echo "${R_SAMPLE_IDS_ARR[*]}")
R_PASSED_SAMPLE_CONDITIONS=$(IFS=,; echo "${R_SAMPLE_CONDITIONS_ARR[*]}")
R_PASSED_SAMPLE_REPLICATES=$(IFS=,; echo "${R_SAMPLE_REPLICATES_ARR[*]}")

# Define paths to consensus peak files
GFP_CONSENSUS_PEAK="${CONSENSUS_PEAKS_DIR}/GFP_consensus_peaks.bed"
YAF_CONSENSUS_PEAK="${CONSENSUS_PEAKS_DIR}/YAF_consensus_peaks.bed"

# Validate all input files first
log_message "Validating input files..."

# Validate consensus peak files
if ! validate_consensus_peak_file "$GFP_CONSENSUS_PEAK"; then
    log_message "ERROR: GFP consensus peak file validation failed."
    exit 1
fi
if ! validate_consensus_peak_file "$YAF_CONSENSUS_PEAK"; then
    log_message "ERROR: YAF consensus peak file validation failed."
    exit 1
fi

# Validate BAM files for selected samples
for sample_id in "${samples_to_process[@]}"; do
    bam_file="${ALIGN_DIR}/${sample_id}.dedup.bam"
    if [[ ! -f "$bam_file" ]]; then
        log_message "ERROR: BAM file not found: $bam_file for sample $sample_id"
        exit 1
    fi
    # Check for BAM index
    if [[ ! -f "${bam_file}.bai" ]]; then
        log_message "WARNING: BAM index file not found for $bam_file. The R script will attempt to create it."
    fi
done
log_message "All required input files appear to be present."

# Run R script for differential binding analysis
log_message "Running differential binding analysis using consensus peaks..."
R_SCRIPT_PATH="${SCRIPT_DIR}/7_differential_binding_consensus_v2.R"

if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    log_message "ERROR: R script not found at $R_SCRIPT_PATH"
    exit 1
fi

# Redirect stdout and stderr to log files
if ! Rscript "$R_SCRIPT_PATH" \
    "${OUTPUT_DIR}" \
    "${ALIGN_DIR}" \
    "${R_PASSED_SAMPLE_IDS}" \
    "${R_PASSED_SAMPLE_CONDITIONS}" \
    "${R_PASSED_SAMPLE_REPLICATES}" \
    "${GFP_CONSENSUS_PEAK}" \
    "${YAF_CONSENSUS_PEAK}" > "${LOG_DIR}/7_differential_binding_consensus_v2.out" 2> "${LOG_DIR}/7_differential_binding_consensus_v2.err"; then
    log_message "ERROR: Differential binding analysis failed. Check logs ${LOG_DIR}/7_differential_binding_consensus_v2.err and .out"
    exit 1
fi

log_message "Differential binding analysis with consensus peaks completed successfully"
log_message "Output files are in ${OUTPUT_DIR}"
log_message "Logs are in ${LOG_DIR}"