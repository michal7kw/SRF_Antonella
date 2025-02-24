#!/bin/bash
#SBATCH --job-name=12_go_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/12_go_enrichment.err"
#SBATCH --output="logs/12_go_enrichment.out"

# Documentation:
# This script performs GO enrichment analysis for different gene sets derived from
# overlap analysis of ChIP-seq data. It processes gene lists to identify enriched
# GO terms, providing insights into the biological functions associated with the genes.
#
# Input files:
# - analysis/11_gene_overlap_analysis/{gene_set}.txt: Text files containing gene lists for GO enrichment.
#
# Output files:
# - analysis/12_go_enrichment/go_enrichment_broad/{gene_set}/go_enrichment_results.csv: CSV files containing GO enrichment results for each gene set.
# - analysis/12_go_enrichment/go_enrichment_broad/{gene_set}/go_enrichment_plots.pdf: PDF files containing GO enrichment plots.
#
# Gene sets:
# - all_overlapping_genes: Genes overlapping between YAF and SOX2.
# - all_yaf_only_genes: Genes uniquely associated with YAF.
# - regulatory_overlapping_genes: Regulatory genes overlapping between YAF and SOX2.
# - regulatory_yaf_only_genes: Regulatory genes uniquely associated with YAF.

# Error handling
set -e
set -u
set -o pipefail

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Define constants
PEAK_TYPE="broad"

# Define the main output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/12_go_enrichment"

# Create the main output directory
log_message "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Define output base directory for GO analysis
output_base="${OUTPUT_DIR}/go_enrichment_${PEAK_TYPE}"
log_message "Creating output base directory: $output_base"
mkdir -p "$output_base"

# Input gene lists directory
input_dir="analysis/11_gene_overlap_analysis"

# List of gene sets to analyze
gene_sets=(
    "all_overlapping_genes"
    "all_yaf_only_genes"
    "regulatory_overlapping_genes"
    "regulatory_yaf_only_genes"
)

# Create subdirectories for each analysis
for gene_set in "${gene_sets[@]}"; do
    log_message "Creating subdirectory: $output_base/$gene_set"
    mkdir -p "$output_base/$gene_set"
done

# Run GO enrichment for each gene set
log_message "Starting GO enrichment analysis..."

for gene_set in "${gene_sets[@]}"; do
    log_message "Processing gene set: $gene_set"
    input_file="$input_dir/${gene_set}.txt"
    output_dir="$output_base/$gene_set"

    if [ -f "$input_file" ]; then
        log_message "Input file found: $input_file"
        python scripts/10_go_enrichment.py \
            --input "$input_file" \
            --output-dir "$output_dir" \
            --description "${gene_set}"
    else
        log_message "Warning: Input file $input_file not found. Skipping $gene_set."
    fi
done

log_message "GO enrichment analysis complete. Results are in $output_base/"
