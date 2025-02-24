#!/bin/bash
#SBATCH --job-name=8_advanced_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/8_advanced_analysis.err"
#SBATCH --output="logs/8_advanced_analysis.out"

# Documentation:
# This script performs analysis of ChIP-seq peaks including:
# - Peak width distribution analysis
# - Signal intensity correlation between samples
# - Genomic distribution of peaks relative to genes
# - Signal profile analysis around TSS regions
# - Motif enrichment analysis (for narrow peaks)
# - Peak clustering based on signal patterns
#
# Input files:
# - analysis/7_differential_binding/diffbind_{peak_type}/significant_peaks.rds: GRanges object with differential peaks
# - analysis/8_annotation_and_enrichment/annotation_{peak_type}/peak_annotation.rds: ChIPseeker annotation object
#
# Output files in analysis/10_additional_analysis/:
#   plots/
#     - peak_width_distribution.pdf: Distribution of peak widths
#     - signal_correlation_heatmap.pdf: Correlation between sample signals
#     - genomic_distribution.pdf: Peak distribution relative to genomic features
#     - tss_profile.pdf: Average signal profile around TSS
#     - motif_enrichment.pdf: Enriched sequence motifs (narrow peaks only)
#     - peak_clusters.pdf: Clustering of peaks by signal patterns
#   summary_statistics.txt: Key metrics from the analysis
#
# Dependencies:
# - GenomicRanges for genomic interval operations
# - ComplexHeatmap and circlize for heatmap visualization
# - ggplot2 for plotting
# - ChIPseeker for genomic feature annotation
# - motifmatchr and JASPAR2020 for motif analysis
# - BSgenome.Hsapiens.UCSC.hg38 for genome sequence
# - DiffBind for peak analysis

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to print usage
usage() {
    echo "Usage: $0 [-t <peak_type>]"
    echo "Options:"
    echo "  -t    Peak type (narrow or broad, default: broad)"
    exit 1
}

# Default peak type
PEAK_TYPE="broad"

# Parse command line options
while getopts "t:h" opt; do
    case $opt in
        t)
            PEAK_TYPE=$OPTARG
            if [[ ! "$PEAK_TYPE" =~ ^(narrow|broad)$ ]]; then
                log_message "ERROR: Peak type must be either 'narrow' or 'broad'"
                usage
            fi
            ;;
        h)
            usage
            ;;
        \?)
            log_message "ERROR: Invalid option -$OPTARG"
            usage
            ;;
    esac
done

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/10_additional_analysis"

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs
mkdir -p ${OUTPUT_DIR}/plots

# Check for required input files
log_message "Checking input files..."
required_files=(
    "analysis/7_differential_binding/diffbind_${PEAK_TYPE}/significant_peaks.rds"
    "analysis/8_annotation_and_enrichment/annotation_${PEAK_TYPE}/peak_annotation.rds"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Run R script for advanced analysis
log_message "Running advanced analysis for ${PEAK_TYPE} peaks..."
Rscript scripts/8_advanced_analysis.R "$PEAK_TYPE" "${OUTPUT_DIR}"

log_message "Advanced analysis completed"
