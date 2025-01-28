#!/bin/bash
#SBATCH --job-name=cross_analysis
#SBATCH --output=logs/cross_analysis.out
#SBATCH --error=logs/cross_analysis.err
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
mkdir -p logs results/{plots,tables,downstream_analysis}

# Activate conda environment
log_message "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/cross_analysis"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Check for required input files
log_message "Checking required input files..."

# V5 peaks and bigwig files
V5_NARROW_PEAKS="../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_narrow_peaks.narrowPeak"
V5_BROAD_PEAKS="../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_broad_peaks.broadPeak"
V5_BIGWIG="../SRF_V5/bigwig/SES-V5ChIP-Seq2_S6.bw"
V5_INPUT_BIGWIG="../SRF_V5/bigwig/InputSES-V5ChIP-Seq_S2.bw"

# H2AK119Ub analysis files and bigwig files
H2A_BASE="../SRF_H2AK119Ub/1_iterative_processing/analysis"
H2A_PEAKS_DIR="${H2A_BASE}/peaks"
H2A_BIGWIG_DIR="${H2A_BASE}/visualization"

required_files=(
    "$V5_NARROW_PEAKS"
    "$V5_BROAD_PEAKS"
    "$V5_BIGWIG"
    "$V5_INPUT_BIGWIG"
    "${H2A_PEAKS_DIR}/GFP_1_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/GFP_2_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/GFP_3_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/YAF_1_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/YAF_2_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/YAF_3_narrow_peaks.narrowPeak"
    "${H2A_PEAKS_DIR}/GFP_1_broad_peaks.broadPeak"
    "${H2A_PEAKS_DIR}/GFP_2_broad_peaks.broadPeak"
    "${H2A_PEAKS_DIR}/GFP_3_broad_peaks.broadPeak"
    "${H2A_PEAKS_DIR}/YAF_1_broad_peaks.broadPeak"
    "${H2A_PEAKS_DIR}/YAF_2_broad_peaks.broadPeak"
    "${H2A_PEAKS_DIR}/YAF_3_broad_peaks.broadPeak"
    "${H2A_BIGWIG_DIR}/GFP_1.bw"
    "${H2A_BIGWIG_DIR}/GFP_2.bw"
    "${H2A_BIGWIG_DIR}/GFP_3.bw"
    "${H2A_BIGWIG_DIR}/YAF_1.bw"
    "${H2A_BIGWIG_DIR}/YAF_2.bw"
    "${H2A_BIGWIG_DIR}/YAF_3.bw"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
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
    "BSgenome.Hsapiens.UCSC.hg38",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "ChIPseeker",
    "clusterProfiler",
    "ComplexHeatmap"
), update = FALSE, ask = FALSE)
EOF

# Run the initial cross-reference analysis
log_message "Starting cross-reference analysis..."
Rscript cross_reference_analysis.R

# Check if the initial analysis was successful
if [ $? -eq 0 ] && [ -f "results/processed_peaks.rds" ]; then
    # Run the downstream analysis
    log_message "Starting downstream analysis..."
    Rscript downstream_analysis.R
else
    log_message "ERROR: Initial cross-reference analysis failed or did not produce required output"
    exit 1
fi

# Check output files from both analyses
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

# Check main results directory
check_directory_contents "results" ""

# Check downstream analysis directory
check_directory_contents "results/downstream_analysis" "downstream/"

log_message "Analysis pipeline completed successfully"
