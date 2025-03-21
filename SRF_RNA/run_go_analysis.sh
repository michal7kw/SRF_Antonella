#!/bin/bash
#SBATCH --job-name=GO_analysis
#SBATCH --output=logs/GO_analysis.out
#SBATCH --error=logs/GO_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Exit on error
set -e

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Run the GO analysis script
python go_analysis.py

echo "Analysis complete! Results are in the results/go_analysis directory." 