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

# Create necessary directories
log_message "Creating output directories..."
mkdir -p logs
mkdir -p analysis/advanced_analysis_broad/{motifs,clusters,profiles,plots}
# mkdir -p analysis/advanced_analysis_narrow/{motifs,clusters,profiles,plots}


# Check for required input files
log_message "Checking input files..."
required_files=(
    "analysis/diffbind_${PEAK_TYPE}/significant_peaks.rds"
    "analysis/annotation_${PEAK_TYPE}/peak_annotation.rds"
)

for file in "${required_files[@]}"; do
    if [[ ! -f "$file" ]]; then
        log_message "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Run R script for advanced analysis
log_message "Running advanced analysis for ${PEAK_TYPE} peaks..."
Rscript scripts/8_advanced_analysis.R "$PEAK_TYPE"

log_message "Advanced analysis completed"
