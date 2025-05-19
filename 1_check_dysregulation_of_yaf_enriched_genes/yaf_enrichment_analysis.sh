#!/bin/bash
#SBATCH --job-name=YAF_binding
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/YAF_binding.err"
#SBATCH --output="logs/YAF_binding.out"

# Load necessary modules (adjust as needed for your cluster)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake 

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5

# Create logs directory if it doesn't exist
mkdir -p logs

Rscript yaf_enrichment_analysis.R