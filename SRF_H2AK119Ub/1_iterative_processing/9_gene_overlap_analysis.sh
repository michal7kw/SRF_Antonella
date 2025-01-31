#!/bin/bash
#SBATCH --job-name=9_gene_overlap_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/9_gene_overlap_analysis.err"
#SBATCH --output="logs/9_gene_overlap_analysis.out"

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
mkdir -p analysis/overlap_analysis_broad logs
# mkdir -p analysis/overlap_analysis_narrow logs

# Run python script for gene overlap analysis
log_message "Running gene overlap analysis..."
python scripts/9_gene_overlap_analysis.py

log_message "Gene overlap analysis completed" 