#!/bin/bash
#SBATCH --job-name=9_advanced_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/9_advanced_analysis.err"
#SBATCH --output="logs/9_advanced_analysis.out"

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
mkdir -p analysis/advanced_analysis/{motifs,clusters,profiles} logs

# Check for required input files
log_message "Checking input files..."
required_files=(
    "analysis/diffbind_broad/significant_peaks.rds"
    "analysis/diffbind_narrow/significant_peaks.rds"
    "analysis/annotation_broad/peak_annotation.rds"
    "analysis/annotation_narrow/peak_annotation.rds"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Run R script for advanced analysis
log_message "Running advanced analysis..."
Rscript scripts/10_advanced_analysis.R

log_message "Advanced analysis completed"
