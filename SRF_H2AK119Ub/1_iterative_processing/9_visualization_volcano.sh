#!/bin/bash
# Documentation:
# This script generates a volcano plot for ChIP-seq peak analysis results.
# It is a simplified version focused solely on the volcano plot.
#
# Input files:
# - A CSV file containing promoter-associated gene annotations.
#   (e.g., analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_promoters.csv)
#
# Output files in analysis/9_visualization_volcano/{PEAK_TYPE}_promoters/peak_analysis/:
#   - volcano_plot_promoters.pdf: Volcano plot of significance vs fold change
#
# Dependencies:
# - Rscript
# - The R script 'scripts/9_visualization_volcano.R'
# - An R environment with ggplot2 and ggrepel packages.
# - 'scripts/utils.R' (expected by the R script)

set -e
set -u
set -o pipefail

# --- Configuration ---
# Define the input CSV file containing promoter gene data
# This file should have 'fold_change', 'FDR', and 'SYMBOL' columns.
INPUT_CSV="analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_promoters.csv"

# Define the peak type (e.g., "broad", "narrow"). This is used for naming output paths.
PEAK_TYPE="broad"

# Define the base output directory for visualizations
OUTPUT_DIR_BASE="analysis/9_visualization_volcano"

# Define the R script to be executed
R_SCRIPT_PATH="scripts/9_visualization_volcano.R"
# --- End Configuration ---

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Ensure the correct conda environment is active if needed.
# Example: conda activate your_r_env
# log_message "Ensure your R environment (e.g., conda) with ggplot2 and ggrepel is active."

# Define working directory (assuming script is run from project root)
WORKDIR="."
cd "$WORKDIR" || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }
log_message "Working directory: $(pwd)"

# Define the specific output directory for this run
# The R script will create subdirectories like {PEAK_TYPE}_promoters/peak_analysis within this.
TARGET_OUTPUT_DIR="${OUTPUT_DIR_BASE}" # R script handles subfolder creation based on peak_type

# Create logs directory and the base output directory if they don't exist
mkdir -p logs "$TARGET_OUTPUT_DIR"
log_message "Base output directory set to: $TARGET_OUTPUT_DIR"

# Check for required R script
if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    log_message "ERROR: R script not found: $R_SCRIPT_PATH"
    exit 1
fi
log_message "Using R script: $R_SCRIPT_PATH"

# Check for required input CSV file
log_message "Checking for input CSV file: $INPUT_CSV"
if [[ ! -f "$INPUT_CSV" ]]; then
    log_message "ERROR: Required input CSV file not found: $INPUT_CSV"
    exit 1
fi
log_message "Input CSV file found: $INPUT_CSV"

# Run R script for volcano plot visualization
log_message "Running volcano plot generation..."
Rscript "${R_SCRIPT_PATH}" "${TARGET_OUTPUT_DIR}" "${INPUT_CSV}" "${PEAK_TYPE}"

log_message "Volcano plot generation completed. Check output in ${TARGET_OUTPUT_DIR}/${PEAK_TYPE}_promoters/peak_analysis/"