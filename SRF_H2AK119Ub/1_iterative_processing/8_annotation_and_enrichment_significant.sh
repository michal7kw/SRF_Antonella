#!/bin/bash

# Documentation:
# This version is adapted for local execution.
# This script performs annotation and enrichment analysis on peaks identified from
# differential binding analysis. It takes the peaks, annotates them relative to genes,
# performs GO enrichment analysis, and generates various visualizations.
#
# Input files:
# - analysis/7_differential_binding_v2/differential_peaks_fdr_lfc_significant.csv: CSV file with significant differential peaks
#
# Output files:
# In analysis/8_annotation_and_enrichment/annotation/:
#   figures/
#     - annotation_plots.pdf: Peak annotation visualizations (pie chart, TSS distance)
#     - detailed_pie_chart.pdf: Detailed pie chart with genomic feature distribution
#     - go_enrichment_plots.pdf: GO term enrichment plots (dotplot, emap, cnet)
#   tables/
#     - peak_annotation.csv: Detailed peak annotations
#     - go_enrichment.csv: GO enrichment analysis results
#   peak_annotation.rds: R object with full annotation data
#
# In analysis/8_annotation_and_enrichment/gene_lists/:
#   - YAF_enriched_genes_full.csv: All enriched genes with details
#   - YAF_enriched_genes_symbols.txt: List of gene symbols only
#   - YAF_enriched_genes_promoters.txt: Genes associated with promoter regions
#   - YAF_enriched_genes_promoters_1st_exon_intron.txt: Genes in promoters + 1st exon/intron
#   - YAF_enriched_genes_distal_intergenic.txt: Genes associated with distal intergenic regions
#   - YAF_enriched_genes_other_regions.txt: Genes in other genomic regions
#
# Dependencies:
# - ChIPseeker for peak annotation
# - clusterProfiler for GO enrichment
# - org.Hs.eg.db for gene ID mapping
# - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
# - ggupset (optional) for upset plots
# - ggplot2 for custom visualizations

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Define working directory
WORKDIR="."
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define directories and input file
INPUT_CSV_FILE="${WORKDIR}/analysis/7_differential_binding_v2/differential_peaks_fdr_lfc_significant.csv"
OUTPUT_DIR="${WORKDIR}/analysis/8_annotation_and_enrichment_significant"

# INPUT_CSV_FILE="${WORKDIR}/analysis/7_differential_binding_strict_v2/differential_peaks_fdr_lfc_significant.csv"
# OUTPUT_DIR="${WORKDIR}/analysis/8_annotation_and_enrichment_strict_significant"

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs

# Check for required input files
log_message "Checking input files..."
if [[ ! -f "$INPUT_CSV_FILE" ]]; then
    log_message "ERROR: Required input CSV file not found: $INPUT_CSV_FILE"
    exit 1
fi

# Run R script for annotation and enrichment analysis
log_message "Running annotation and enrichment analysis..."
Rscript scripts/8_annotation_and_enrichment_significant.R "${OUTPUT_DIR}" "${INPUT_CSV_FILE}"

log_message "Annotation and enrichment analysis completed"

# dos2unix SRF_H2AK119Ub/1_iterative_processing/8_annotation_and_enrichment_significant.sh

