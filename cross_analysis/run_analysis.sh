#!/bin/bash
#SBATCH --job-name=cross_ref
#SBATCH --output=logs/cross_ref.out
#SBATCH --error=logs/cross_ref.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/cross_analysis

Rscript cross_reference_analysis.R
