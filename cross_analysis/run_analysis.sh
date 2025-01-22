#!/bin/bash
#SBATCH --job-name=cross_ref
#SBATCH --output=logs/cross_ref.out
#SBATCH --error=logs/cross_ref.err
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
mkdir -p logs results/{plots,tables}

# Activate conda environment
log_message "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/cross_analysis"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Check for required input files
log_message "Checking required input files..."

# V5 peaks
V5_PEAKS="../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_peaks.narrowPeak"
if [[ ! -f "$V5_PEAKS" ]]; then
    log_message "ERROR: V5 peaks file not found: $V5_PEAKS"
    exit 1
fi

# H2AK119Ub analysis files
H2A_BASE="../SRF_H2AK119Ub/1_iterative_processing/analysis"
required_files=(
    "${H2A_BASE}/diffbind_narrow/differential_peaks.csv"
    "${H2A_BASE}/diffbind_broad/differential_peaks.csv"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
    log_message "Found input file: $file"
done

# Run the analysis
log_message "Starting cross-reference analysis..."
Rscript cross_reference_analysis.R

# Check output files
log_message "Checking output files..."
required_outputs=(
    "results/tables/overlap_statistics.csv"
    "results/tables/overlapping_genes.csv"
    "results/plots/overlap_comparison.pdf"
    "results/plots/venn_diagram_narrow.png"
    "results/plots/venn_diagram_broad.png"
    "results/plots/gene_overlap_venn.png"
    "results/analysis_summary.txt"
)

for file in "${required_outputs[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "WARNING: Expected output file not found: $file"
    else
        log_message "Generated output file: $file ($(du -h "$file" | cut -f1))"
    fi
done

log_message "Analysis completed successfully"
