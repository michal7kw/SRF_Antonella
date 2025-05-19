#!/bin/bash

# Documentation:
# This script filters the 'differential_peaks.csv' file (output from DiffBind analysis)
# to retain statistically significant peaks.
# It first filters by a specified FDR threshold.
# Optionally, it can then further filter those peaks by a log fold change (LFC) threshold,
# saving the results to two separate files.
# It uses the R script '7b_filter_significant_peaks.R'.

# Error handling
set -e  # Exit immediately if a command exits with non-zero status
set -u  # Treat unset variables as an error
set -o pipefail  # Fail pipeline if any command fails

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Define working directory (assuming script is run from 1_iterative_processing)
WORKDIR="."
cd "$WORKDIR" || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define directories and files
INPUT_CSV_DIR="analysis/7_differential_binding_v2"
# INPUT_CSV_DIR="analysis/7_differential_binding_strict_v2"

INPUT_CSV_FILENAME="differential_peaks.csv"
INPUT_CSV_PATH="${INPUT_CSV_DIR}/${INPUT_CSV_FILENAME}"

# Output directory and CSV filename for FDR significant peaks
OUTPUT_FDR_DIR="${INPUT_CSV_DIR}" # Or e.g., "analysis/7b_filtered_peaks"
OUTPUT_FDR_FILENAME="differential_peaks_fdr_significant.csv"
OUTPUT_FDR_CSV_PATH="${OUTPUT_FDR_DIR}/${OUTPUT_FDR_FILENAME}"

# Output directory and CSV filename for FDR AND LFC significant peaks
# This is the new second output file.
OUTPUT_LFC_DIR="${INPUT_CSV_DIR}" # Can be the same or different directory
OUTPUT_LFC_FILENAME="differential_peaks_fdr_lfc_significant.csv"
OUTPUT_LFC_CSV_PATH="${OUTPUT_LFC_DIR}/${OUTPUT_LFC_FILENAME}"


# R script for filtering
R_SCRIPT_PATH="scripts/7b_filter_significant_peaks.R"

# Parameters for the R script
FDR_THRESHOLD=0.05  # Common threshold for FDR significance
LFC_THRESHOLD=1.0   # Absolute log2 fold change threshold (e.g., 1.0 for a 2-fold change)
LFC_COLUMN_NAME="Fold" # Column name in the CSV that contains LFC values

# Control whether to generate the LFC filtered file. Set to "true" to enable.
# If "false", only the FDR filtered file will be generated (or attempted).
GENERATE_LFC_FILE="true"


# Create output directory if it doesn't exist (if different from input)
# mkdir -p "${OUTPUT_FDR_DIR}" # Uncomment if OUTPUT_FDR_DIR is new
# mkdir -p "${OUTPUT_LFC_DIR}" # Uncomment if OUTPUT_LFC_DIR is new

log_message "Starting filtering of differential peaks."

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    log_message "ERROR: Rscript command not found. Please ensure R is installed and in your PATH."
    exit 1
fi

# Check if the R filtering script exists
if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    log_message "ERROR: R filtering script not found at $R_SCRIPT_PATH"
    exit 1
fi

# Check if the input CSV file exists
if [[ ! -f "$INPUT_CSV_PATH" ]]; then
    log_message "ERROR: Input CSV file not found at $INPUT_CSV_PATH"
    log_message "Please ensure the DiffBind analysis script has run successfully to generate ${INPUT_CSV_FILENAME}."
    exit 1
fi

log_message "Input CSV (all peaks): ${INPUT_CSV_PATH}"
log_message "Output CSV (FDR significant peaks): ${OUTPUT_FDR_CSV_PATH}"
log_message "FDR Threshold: ${FDR_THRESHOLD}"

if [[ "$GENERATE_LFC_FILE" == "true" ]]; then
    log_message "Output CSV (FDR & LFC significant peaks): ${OUTPUT_LFC_CSV_PATH}"
    log_message "LFC Threshold (absolute): ${LFC_THRESHOLD} (using column '${LFC_COLUMN_NAME}')"
fi
log_message "R script: ${R_SCRIPT_PATH}"


# Construct Rscript command
R_COMMAND="Rscript ${R_SCRIPT_PATH} \
    --input_csv \"${INPUT_CSV_PATH}\" \
    --output_csv_fdr \"${OUTPUT_FDR_CSV_PATH}\" \
    --fdr_threshold ${FDR_THRESHOLD} \
    --lfc_column_name \"${LFC_COLUMN_NAME}\""

if [[ "$GENERATE_LFC_FILE" == "true" ]]; then
    R_COMMAND="${R_COMMAND} --output_csv_lfc \"${OUTPUT_LFC_CSV_PATH}\" --lfc_threshold ${LFC_THRESHOLD}"
fi

# Execute the R script
log_message "Running R script for filtering peaks..."
if ! eval "$R_COMMAND"; then
    log_message "ERROR: R script execution failed. Check R script output for details."
    exit 1
fi

log_message "Filtering for peaks completed successfully."
log_message "Output CSV file with FDR significant peaks: ${OUTPUT_FDR_CSV_PATH}"
if [[ "$GENERATE_LFC_FILE" == "true" ]]; then
    log_message "Output CSV file with FDR & LFC significant peaks: ${OUTPUT_LFC_CSV_PATH}"
    log_message "The second CSV can be used for analyses requiring stricter criteria."
fi
log_message ""
log_message "These CSV files can now be used as input for downstream analyses."

log_message "Script finished."