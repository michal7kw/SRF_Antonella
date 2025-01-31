#!/bin/bash
#SBATCH --job-name=5_differential_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding.err"
#SBATCH --output="logs/5_differential_binding.out"

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
mkdir -p analysis/diffbind_broad
# mkdir -p analysis/diffbind_narrow
mkdir -p analysis/plots_broad
# mkdir -p analysis/plots_narrow


# Check for required input files
log_message "Checking input files..."
for sample in GFP_{1..3} YAF_{1..3}; do
    for type in broad; do #in broad narrow; do
        files=(
            "analysis/aligned/${sample}.dedup.bam"
            "analysis/peaks/${sample}_${type}_peaks.${type}Peak"
        )
        for file in "${files[@]}"; do
            if [[ ! -f "$file" ]]; then
                log_message "ERROR: Required file not found: $file"
                exit 1
            fi
        done
    done
done

# Run R script for differential binding analysis
log_message "Running differential binding analysis..."
Rscript scripts/5_differential_binding.R

log_message "Differential binding analysis completed" 