#!/bin/bash
#SBATCH --job-name=5_differential_binding_broad
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_differential_binding_broad.err"
#SBATCH --output="logs/5_differential_binding_broad.out"

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check if output exists and is not empty
check_output() {
    local file=$1
    if [[ ! -s $file ]]; then
        log_message "ERROR: Output file $file is empty or does not exist"
        exit 1
    else
        log_message "Successfully created: $file ($(du -h $file | cut -f1))"
    fi
}

# Function to check if a command exists
check_command() {
    if ! command -v $1 &> /dev/null; then
        log_message "ERROR: $1 command not found. Please ensure it is installed and in PATH"
        exit 1
    fi
}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Check required commands
log_message "Checking required commands..."
check_command R
check_command Rscript

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/{diffbind_broad,annotation_broad,qc} logs || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check for required input files
log_message "Checking required input files..."
required_files=(
    "analysis/peaks/GFP_1_broad_peaks.broadPeak"
    "analysis/peaks/GFP_2_broad_peaks.broadPeak"
    "analysis/peaks/GFP_3_broad_peaks.broadPeak"
    "analysis/peaks/YAF_1_broad_peaks.broadPeak"
    "analysis/peaks/YAF_2_broad_peaks.broadPeak"
    "analysis/peaks/YAF_3_broad_peaks.broadPeak"
    "analysis/aligned/GFP_1.dedup.bam"
    "analysis/aligned/GFP_2.dedup.bam"
    "analysis/aligned/GFP_3.dedup.bam"
    "analysis/aligned/YAF_1.dedup.bam"
    "analysis/aligned/YAF_2.dedup.bam"
    "analysis/aligned/YAF_3.dedup.bam"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
    log_message "Found input file: $file ($(du -h $file | cut -f1))"
done

# Check for required R packages
log_message "Checking R package dependencies..."
R --quiet -e '
    packages <- c("DiffBind", "ChIPseeker", "clusterProfiler", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
    missing <- packages[!packages %in% installed.packages()[,"Package"]]
    if (length(missing) > 0) {
        message("Installing missing packages: ", paste(missing, collapse=", "))
        if (!require("BiocManager", quietly=TRUE)) {
            install.packages("BiocManager", repos="https://cloud.r-project.org")
        }
        BiocManager::install(missing, update=FALSE)
    }
    sapply(packages, require, character.only=TRUE)
' || { log_message "ERROR: Failed to install/load required R packages"; exit 1; }

# Run differential binding analysis
log_message "Starting differential binding analysis for broad peaks..."
Rscript scripts/5_differential_binding_broad.R

# Check output files
log_message "Checking output files..."
output_files=(
    "analysis/diffbind_broad/all_peaks.rds"
    "analysis/diffbind_broad/significant_peaks.rds"
    "analysis/diffbind_broad/differential_peaks.csv"
    "analysis/annotation_broad/peak_annotation.rds"
    "analysis/annotation_broad/peak_annotation.csv"
    "analysis/annotation_broad/annotation_plots.pdf"
)

for file in "${output_files[@]}"; do
    check_output "$file"
done

log_message "Differential binding analysis for broad peaks completed successfully"
