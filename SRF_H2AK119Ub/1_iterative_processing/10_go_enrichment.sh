#!/bin/bash

# Stop on error
set -e

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Create output directory for GO analysis
output_base="analysis/go_enrichment"
mkdir -p "$output_base"

# Input gene lists directory
input_dir="analysis/overlap_analysis"

# List of gene sets to analyze
gene_sets=(
    "all_overlapping_genes"
    "all_yaf_only_genes"
    "regulatory_overlapping_genes"
    "regulatory_yaf_only_genes"
)

# Create subdirectories for each analysis
for gene_set in "${gene_sets[@]}"; do
    mkdir -p "$output_base/$gene_set"
done

# Run GO enrichment for each gene set
log_message "Starting GO enrichment analysis..."

for gene_set in "${gene_sets[@]}"; do
    log_message "Processing $gene_set..."
    input_file="$input_dir/${gene_set}.txt"
    output_dir="$output_base/$gene_set"
    
    if [ -f "$input_file" ]; then
        python scripts/10_go_enrichment.py \
            --input "$input_file" \
            --output-dir "$output_dir" \
            --description "${gene_set}"
    else
        log_message "Warning: Input file $input_file not found"
    fi
done

log_message "GO enrichment analysis complete. Results are in $output_base/"
