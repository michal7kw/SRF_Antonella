#!/bin/bash
#SBATCH --job-name=5_differential_binding_parallel
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding_parallel.err"
#SBATCH --output="logs/5_differential_binding_parallel.out"
#SBATCH --array=0-5   # Added array job for 6 samples

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

PEAKS_SUFFIX="peaks_final"
# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs
mkdir -p analysis/diffbind_broad_parallel
mkdir -p analysis/plots_broad_parallel

# Define sample list and select the current sample
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample="${samples[$SLURM_ARRAY_TASK_ID]}"
log_message "Processing sample: ${sample}"

# For this sample, check for the required input files (only using type 'broad')
type="broad"
files=(
    "analysis/aligned/${sample}.dedup.bam"
    "analysis/peaks/${sample}_${type}_${PEAKS_SUFFIX}.${type}Peak"
)
log_message "Checking input files for ${sample}..."
for file in "${files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Run R script for differential binding analysis for the current sample
log_message "Running differential binding analysis for ${sample}..."
Rscript scripts/5_differential_binding.R "${sample}"

log_message "Differential binding analysis completed for ${sample}" 