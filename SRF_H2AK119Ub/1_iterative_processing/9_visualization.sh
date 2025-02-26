#!/bin/bash
#SBATCH --job-name=9_visualization
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/9_visualization.err"
#SBATCH --output="logs/9_visualization.out"

# Documentation:
# This script creates visualizations for ChIP-seq peak analysis results.
# It generates plots and statistics to help interpret the differential binding analysis
# between YAF and GFP samples.
#
# Input files:
# - analysis/7_differential_binding/diffbind_broad/significant_peaks.rds: GRanges object with differential peaks
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_full.csv: Annotated gene list
# - analysis/7_differential_binding/diffbind_narrow/significant_peaks.rds (optional): Narrow peaks data
# - analysis/8_annotation_and_enrichment/gene_lists_narrow/YAF_enriched_genes_narrow_full.csv (optional): Narrow peaks genes
#
# Output files in analysis/9_visualization:
#   peak_analysis/
#     - peak_distribution.pdf: Bar plot showing distribution of peak categories
#     - ma_plot_enhanced.pdf: MA plot with fold change vs concentration
#     - volcano_plot_enhanced.pdf: Volcano plot of significance vs fold change
#   gene_analysis/
#     - tss_distribution.pdf: Histogram of peak distances to TSS
#   summary_statistics/
#     - summary_stats.csv: Table with key analysis metrics
#
# Dependencies:
# - ggplot2 for plotting
# - gridExtra for plot arrangement
# - scales for axis formatting
# - utils.R for helper functions

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
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/9_visualization"
mkdir -p logs

# Check for required input files
log_message "Checking input files..."
for type in broad; do
    files=(
        "analysis/7_differential_binding/significant_peaks.rds"
        "analysis/8_annotation_and_enrichment/gene_lists_${type}/YAF_enriched_genes_${type}_full.csv"
    )
    for file in "${files[@]}"; do
        if [[ ! -f "$file" ]]; then
            log_message "ERROR: Required file not found: $file"
            exit 1
        fi
    done
done

# Run R script for visualization
log_message "Running visualization analysis..."
Rscript scripts/9_visualization.R ${OUTPUT_DIR}

log_message "Visualization analysis completed"