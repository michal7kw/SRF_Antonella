#!/bin/bash
#SBATCH --job-name=cross_analysis
#SBATCH --output=logs/cross_analysis.out
#SBATCH --error=logs/cross_analysis.err
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Create necessary directories
log_message "Creating directories..."
mkdir -p logs results/plots

# Activate conda environment
log_message "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/cross_analysis"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Check for required input files
log_message "Checking required input files..."


# Run the analysis
log_message "Starting cross-reference analysis..."
Rscript cross_reference_analysis.R

log_message "Analysis completed successfully"
