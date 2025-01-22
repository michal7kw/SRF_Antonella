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
mkdir -p analysis/{gene_lists_narrow,gene_lists_broad,gene_lists_combined,plots_narrow,plots_broad} logs || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check for required input files for narrow peaks
log_message "Checking required input files for narrow peaks..."
required_files_narrow=(
    "analysis/diffbind_narrow/differential_peaks.csv"
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
    packages <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
                 "clusterProfiler", "dplyr")
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

# Run YAF-enriched genes analysis
log_message "Starting YAF-enriched genes analysis for both narrow and broad peaks..."
Rscript scripts/get_YAF_enriched_genes.R

# Check output files for narrow peaks
log_message "Checking output files for narrow peaks..."
output_files_narrow=(
    "analysis/gene_lists_narrow/YAF_enriched_genes_narrow_full.csv"
    "analysis/gene_lists_narrow/YAF_enriched_genes_narrow_symbols.txt"
    "analysis/gene_lists_narrow/GO_enrichment_narrow.csv"
    "analysis/plots_narrow/peak_annotation_narrow.pdf"
    "analysis/plots_narrow/GO_enrichment_narrow.pdf"
)

# Check output files for broad peaks
log_message "Checking output files for broad peaks..."
output_files_broad=(
    "analysis/gene_lists_broad/YAF_enriched_genes_broad_full.csv"
    "analysis/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt"
    "analysis/gene_lists_broad/GO_enrichment_broad.csv"
    "analysis/plots_broad/peak_annotation_broad.pdf"
    "analysis/plots_broad/GO_enrichment_broad.pdf"
)

# Check combined analysis files
log_message "Checking combined analysis files..."
output_files_combined=(
    "analysis/gene_lists_combined/combined_analysis_summary.txt"
)

# Check all output files
for file in "${output_files_narrow[@]}" "${output_files_broad[@]}" "${output_files_combined[@]}"; do
    check_output "$file"
done

log_message "YAF-enriched genes analysis completed successfully for both narrow and broad peaks"