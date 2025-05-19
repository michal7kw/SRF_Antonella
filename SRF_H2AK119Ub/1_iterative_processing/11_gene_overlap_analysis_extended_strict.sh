#!/bin/bash
# This script runs the extended gene overlap analysis.

#### Input files (handled by the Python script): ####
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_promoters.txt
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_full.csv
# - COMMON_DATA/sox2_binding.csv
# - SRF_RNA/results/deseq2/YAF_vs_GFP/summary_files/YAF_vs_GFP/down_regulated.csv

#### Output files (in ./SRF_H2AK119Ub/1_iterative_processing/analysis/11_gene_overlap_analysis_extended_strict/ relative to project root): ####
# See Python script header for detailed list of outputs.
# Includes 2-way and 3-way Venn diagrams, gene lists, and enrichment stats/plots.

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Define output directory relative to project root
OUTPUT_DIR_PARAM="./analysis/11_gene_overlap_analysis_extended_strict"

# Create necessary directories
log_message "Creating output directories..."
mkdir -p "${OUTPUT_DIR_PARAM}"
# Assuming a 'logs' directory at the project root for general logs, if needed by other parts of a larger workflow.
# If specific logs for this script are desired, they should be handled within the script or Python script.
# mkdir -p logs 

# Run python script for gene overlap analysis
log_message "Running extended gene overlap analysis..."
python ./scripts/11_gene_overlap_analysis_extended_strict.py "${OUTPUT_DIR_PARAM}"

log_message "Extended gene overlap analysis completed. Output in ${OUTPUT_DIR_PARAM}"