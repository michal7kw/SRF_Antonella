#!/bin/bash
#SBATCH --job-name=7_visualization
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_visualization.err"
#SBATCH --output="logs/7_visualization.out"

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

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/plots_{broad,narrow}/{peak_analysis,gene_analysis} logs
mkdir -p analysis/plots_combined

# Check for required input files
log_message "Checking input files..."
for type in broad narrow; do
    files=(
        "analysis/diffbind_${type}/significant_peaks.rds"
        "analysis/gene_lists_${type}/YAF_enriched_genes_${type}_full.csv"
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
Rscript scripts/7_visualization.R

log_message "Visualization analysis completed" 