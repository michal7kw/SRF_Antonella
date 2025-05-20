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
# PLOT_TYPE="Promoters"
# INPUT_CSV="analysis/8_annotation_and_enrichment_strict_all/gene_lists_broad/YAF_enriched_genes_broad_promoters.csv"

PLOT_TYPE="All_Peaks"
INPUT_CSV="analysis/8_annotation_and_enrichment_strict_all/gene_lists_broad/YAF_enriched_genes_broad_full.csv"

# Define the base output directory for visualizations
OUTPUT_DIR_BASE="analysis/9_visualization_volcano_strict"

# Define the R script to be executed
R_SCRIPT_PATH="scripts/9_visualization_volcano.R"
# --- End Configuration ---

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Define working directory
WORKDIR="."
cd "$WORKDIR" || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }
log_message "Working directory: $(pwd)"

# Define the specific output directory for this run
TARGET_OUTPUT_DIR="${OUTPUT_DIR_BASE}"

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
Rscript "${R_SCRIPT_PATH}" "${TARGET_OUTPUT_DIR}" "${INPUT_CSV}" "${PLOT_TYPE}"

log_message "Volcano plot generation completed. Check output in ${TARGET_OUTPUT_DIR}/${PLOT_TYPE}_promoters/peak_analysis/"