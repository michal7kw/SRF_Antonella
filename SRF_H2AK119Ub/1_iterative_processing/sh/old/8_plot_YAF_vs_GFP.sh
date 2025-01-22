#!/bin/bash
#SBATCH --job-name=8_plot_YAF_vs_GFP
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/8_plot_YAF_vs_GFP.err"
#SBATCH --output="logs/8_plot_YAF_vs_GFP.out"

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
mkdir -p analysis/{plots_narrow,plots_broad,plots_combined,gene_lists,qc} logs || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check for required input files for narrow peaks
log_message "Checking required input files for narrow peaks..."
required_files_narrow=(
    "analysis/diffbind_narrow/differential_peaks.csv"
    "analysis/gene_lists/YAF_enriched_genes.txt"
    "analysis/aligned/GFP_1.dedup.bam"
    "analysis/aligned/GFP_2.dedup.bam"
    "analysis/aligned/GFP_3.dedup.bam"
    "analysis/aligned/YAF_1.dedup.bam"
    "analysis/aligned/YAF_2.dedup.bam"
    "analysis/aligned/YAF_3.dedup.bam"
)

# Check for required input files for broad peaks
log_message "Checking required input files for broad peaks..."
required_files_broad=(
    "analysis/diffbind_broad/differential_peaks.csv"
)

# Check all required files
for file in "${required_files_narrow[@]}" "${required_files_broad[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
    log_message "Found input file: $file ($(du -h $file | cut -f1))"
done

# Check for required R packages
log_message "Checking R package dependencies..."
R --quiet -e '
    packages <- c("ggplot2", "dplyr", "tidyr", "GenomicRanges", "rtracklayer", "DiffBind", "ComplexHeatmap", "circlize", "RColorBrewer")
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

# Run visualization script
log_message "Starting visualization for both narrow and broad peaks..."
Rscript scripts/8_plot_YAF_vs_GFP.R

# Check output files for narrow peaks
log_message "Checking output files for narrow peaks..."
output_files_narrow=(
    "analysis/plots_narrow/YAF_vs_GFP_heatmap.pdf"
    "analysis/plots_narrow/YAF_vs_GFP_volcano.pdf"
    "analysis/plots_narrow/YAF_vs_GFP_coverage.pdf"
)

# Check output files for broad peaks
log_message "Checking output files for broad peaks..."
output_files_broad=(
    "analysis/plots_broad/MA_plot_broad.pdf"
    "analysis/plots_broad/volcano_plot_broad.pdf"
    "analysis/plots_broad/summary_stats_broad.txt"
)

# Check combined analysis files
log_message "Checking combined analysis files..."
output_files_combined=(
    "analysis/plots_combined/combined_summary_stats.txt"
)

# Check all output files
for file in "${output_files_narrow[@]}" "${output_files_broad[@]}" "${output_files_combined[@]}"; do
    check_output "$file"
done

log_message "Visualization completed successfully for both narrow and broad peaks"