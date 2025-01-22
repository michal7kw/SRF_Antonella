#!/bin/bash
#SBATCH --job-name=7_get_YAF_enriched_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_get_YAF_enriched_genes.err"
#SBATCH --output="logs/7_get_YAF_enriched_genes.out"

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
mkdir -p analysis/{plots,gene_lists,qc} logs || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check for required input files from previous steps
log_message "Checking required input files..."
required_files=(
    "analysis/diffbind/differential_peaks.csv"
    "analysis/annotation/tables/peak_annotation_full.txt"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        log_message "Please run the differential binding and peak annotation steps first"
        exit 1
    fi
    log_message "Found input file: $file ($(du -h $file | cut -f1))"
done

# Check for required R packages
log_message "Checking R package dependencies..."
R --quiet -e '
    packages <- c("DiffBind", "ChIPseeker", "clusterProfiler", 
                 "org.Hs.eg.db", "ggplot2", "dplyr")
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

# Run YAF enriched genes analysis
log_message "Starting YAF-enriched gene analysis..."
Rscript scripts/get_YAF_enriched_genes.R

# Check output files
log_message "Checking output files..."
output_files=(
    "analysis/gene_lists/YAF_enriched_genes.txt"
    "analysis/gene_lists/YAF_enriched_GO_terms.txt"
    "analysis/plots/YAF_enriched_GO_dotplot.pdf"
)

for file in "${output_files[@]}"; do
    check_output "$file"
done

log_message "YAF-enriched gene analysis completed successfully"