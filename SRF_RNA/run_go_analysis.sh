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
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA || exit 1 # Exit if cd fails

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake || exit 1 # Exit if conda fails

# --- Define parameters ---
DE_FILE="results/deseq2/YAF_vs_GFP/differential_expression.csv"
OUTPUT_DIR="results/go_analysis"
PADJ_THRESHOLD=0.05
FOLD_CHANGE_THRESHOLD=0.0
TOP_N=20
CACHE_FILE="ensembl_symbol_cache.pkl"
# Ensure paths are absolute if needed, or relative to the CWD set above
# Using absolute path for DE_FILE as in original script might be safer
DE_FILE_ABS="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/${DE_FILE}"

# Define libraries as an array
ENRICHR_LIBS=(
    'GO_Biological_Process_2021'
    'GO_Cellular_Component_2021'
    'GO_Molecular_Function_2021'
    'KEGG_2021_Human'
    'Reactome_2022'
)

# Run the GO analysis script with parameters
echo "Starting GO analysis script..."
python go_analysis.py \
  --de_file "${DE_FILE_ABS}" \
  --output_dir "${OUTPUT_DIR}" \
  --padj_threshold ${PADJ_THRESHOLD} \
  --fold_change_threshold ${FOLD_CHANGE_THRESHOLD} \
  --top_n_plots ${TOP_N} \
  --cache_file "${CACHE_FILE}" \
  --enrichr_libs "${ENRICHR_LIBS[@]}" # Pass array elements as separate args

echo "Analysis complete! Results are in the ${OUTPUT_DIR} directory." 