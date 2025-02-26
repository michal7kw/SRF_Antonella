#!/bin/bash
#SBATCH --job-name=8_annotation_and_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/8_annotation_and_enrichment.err"
#SBATCH --output="logs/8_annotation_and_enrichment.out"

# Documentation:
# This script performs annotation and enrichment analysis on peaks identified from
# differential binding analysis. It takes the peaks, annotates them relative to genes,
# performs GO enrichment analysis, and generates various visualizations.
#
# Input files:
# - analysis/7_differential_binding/significant_peaks.rds: GRanges object with significant differential peaks
#
# Output files:
# In analysis/8_annotation_and_enrichment/annotation_broad/:
#   figures/
#     - annotation_plots.pdf: Peak annotation visualizations (pie chart, TSS distance)
#     - detailed_pie_chart.pdf: Detailed pie chart with genomic feature distribution
#     - go_enrichment_plots.pdf: GO term enrichment plots (dotplot, emap, cnet)
#   tables/
#     - peak_annotation.csv: Detailed peak annotations
#     - go_enrichment.csv: GO enrichment analysis results
#   peak_annotation.rds: R object with full annotation data
#
# In analysis/8_annotation_and_enrichment/gene_lists_broad/:
#   - YAF_enriched_genes_broad_full.csv: All enriched genes with details
#   - YAF_enriched_genes_broad_symbols.txt: List of gene symbols only
#   - YAF_enriched_genes_broad_promoters.txt: Genes associated with promoter regions
#   - YAF_enriched_genes_broad_promoters_1st_exon_intron.txt: Genes in promoters + 1st exon/intron
#   - YAF_enriched_genes_broad_distal_intergenic.txt: Genes associated with distal intergenic regions
#   - YAF_enriched_genes_broad_other_regions.txt: Genes in other genomic regions
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

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment"

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs
mkdir -p ${OUTPUT_DIR}/annotation_broad/{figures,tables}
mkdir -p ${OUTPUT_DIR}/gene_lists_broad

# Check for required input files
log_message "Checking input files..."
peak_type="broad"
files=(
    "analysis/7_differential_binding/significant_peaks.rds"
)
for file in "${files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Run R script for annotation and enrichment analysis
log_message "Running annotation and enrichment analysis..."
Rscript scripts/8_annotation_and_enrichment.R ${OUTPUT_DIR}

log_message "Annotation and enrichment analysis completed"