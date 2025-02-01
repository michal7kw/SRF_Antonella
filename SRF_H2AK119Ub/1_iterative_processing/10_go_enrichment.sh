#!/bin/bash
#SBATCH --job-name=10_go_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/10_go_enrichment.err"
#SBATCH --output="logs/10_go_enrichment.out"

# Stop on error
set -e

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

PEAK_TYPE="broad"

# Create output directory for GO analysis
output_base="analysis/go_enrichment_${PEAK_TYPE}"
mkdir -p "$output_base"

# Input gene lists directory
input_dir="analysis/overlap_analysis_${PEAK_TYPE}"

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
