#!/bin/bash
#SBATCH --job-name=6_annotation_and_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/6_annotation_and_enrichment.err"
#SBATCH --output="logs/6_annotation_and_enrichment.out"

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
mkdir -p logs
mkdir -p analysis/annotation_broad/{figures,tables}
# mkdir -p analysis/annotation_narrow/{figures,tables}
mkdir -p analysis/gene_lists_broad
# mkdir -p analysis/gene_lists_narrow

# Check for required input files
log_message "Checking input files..."
for type in broad; do #in broad narrow
    files=(
        "analysis/diffbind_${type}/all_peaks.rds"
        "analysis/diffbind_${type}/significant_peaks.rds"
    )
    for file in "${files[@]}"; do
        if [[ ! -f "$file" ]]; then
            log_message "ERROR: Required file not found: $file"
            exit 1
        fi
    done
done

# Run R script for annotation and enrichment analysis
log_message "Running annotation and enrichment analysis..."
Rscript scripts/6_annotation_and_enrichment.R

log_message "Annotation and enrichment analysis completed" 