#!/bin/bash
#SBATCH --job-name=downstream_analysis
#SBATCH --output=logs/downstream_analysis.out
#SBATCH --error=logs/downstream_analysis.err
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Create necessary directories
log_message "Creating directories..."
mkdir -p logs results/downstream_analysis

# Activate conda environment
log_message "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/cross_analysis"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Check for required input files from previous analysis
log_message "Checking required input files..."
required_inputs=(
    "results/processed_peaks.rds"
    "results/narrow_broad_categorized_peaks.rds"
    "results/broad_broad_categorized_peaks.rds"
    "results/normalized_v5_signal.rds"
)

for file in "${required_inputs[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required input file not found: $file"
        log_message "Please run cross_reference_analysis.R first"
        exit 1
    fi
    log_message "Found input file: $file"
done

# Install required Bioconductor packages
log_message "Installing required Bioconductor packages..."
R --quiet --no-save << EOF
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "ChIPseeker",
    "clusterProfiler",
    "ComplexHeatmap"
), update = FALSE, ask = FALSE)
EOF

# Run the downstream analysis
log_message "Starting downstream analysis..."
Rscript downstream_analysis.R

# Check output files
log_message "Checking output files..."

# Function to check directory contents
check_directory_contents() {
    local dir=$1
    local prefix=$2
    
    if [ -d "$dir" ]; then
        for file in "$dir"/*; do
            if [ -f "$file" ]; then
                log_message "Generated output file: ${prefix}${file#$dir} ($(du -h "$file" | cut -f1))"
            fi
        done
    else
        log_message "WARNING: Directory not found: $dir"
    fi
}

# Check downstream analysis directory
check_directory_contents "results/downstream_analysis" ""

# Check for specific important output files
important_outputs=(
    "results/downstream_analysis/analysis_summary.txt"
    "results/downstream_analysis/broad_broad_v5_with_h2a_genomic_annotation.pdf"
    "results/downstream_analysis/broad_broad_v5_only_genomic_annotation.pdf"
)

for file in "${important_outputs[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "WARNING: Important output file not found: $file"
    fi
done

log_message "Downstream analysis completed successfully" 