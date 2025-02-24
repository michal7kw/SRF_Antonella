#!/bin/bash
#SBATCH --job-name=4_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_multiqc.err"
#SBATCH --output="logs/4_multiqc.out"

# Documentation:
# MultiQC Aggregation and QC Summary Generation
# This script aggregates quality control metrics from multiple pipeline steps using MultiQC
# and generates comprehensive QC reports and summary statistics.

# Key Features:
# - Aggregates FastQC reports from pre- and post-trimming steps
# - Incorporates alignment statistics and duplicate metrics
# - Generates a comprehensive MultiQC report
# - Creates a detailed CSV summary of key QC metrics
# - Handles all 6 samples (GFP_1-3 and YAF_1-3)

# Pipeline Steps:
# 1. Verify required tools and environment
# 2. Collect and aggregate QC data using MultiQC
# 3. Generate detailed QC summary statistics
# 4. Create final reports and metrics files

# Required Tools:
# - MultiQC (v1.14 or later)
# - jq (for JSON parsing)
# - GNU core utilities

# Input Files (per sample):
# Quality Control:
# - analysis/2_quality_control/fastqc/pre_trim/{sample}_R1_001_fastqc.*: Pre-trim FastQC reports
# - analysis/2_quality_control/fastqc/post_trim/{sample}_R1_paired_fastqc.*: Post-trim FastQC reports
# - analysis/2_quality_control/qc/{sample}_trimming_metrics.json: Trimming statistics in JSON format

# Alignment Metrics:
# - analysis/3_alignment/qc/{sample}_alignment_summary.txt: Alignment summary statistics
# - analysis/3_alignment/qc/{sample}_flagstat.txt: SAMtools flagstat output
# - analysis/3_alignment/qc/{sample}_dup_metrics.txt: Picard duplicate metrics
# - analysis/3_alignment/qc/{sample}_insert_size_metrics.txt: Insert size metrics

# Output Files:
# MultiQC Reports:
# - analysis/4_multiqc/multiqc_final.html: Comprehensive MultiQC report
# - analysis/4_multiqc/multiqc_final_data/: MultiQC report data directory
# - logs/multiqc_final.log: MultiQC processing log

# QC Summary:
# - analysis/4_multiqc/qc_summary.csv: CSV file with key QC metrics for all samples
#   Columns: Sample,RawReads,CleanReads,TrimRetentionRate,MappedReads,UniqueReads,FinalRetentionRate

# Exit script on error, undefined variables, and pipe failures
set -e
set -u
set -o pipefail

# Helper function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Helper function to verify output files exist and are not empty
check_output() {
    local file=$1
    if [[ ! -s $file ]]; then
        log_message "ERROR: Output file $file is empty or does not exist"
        exit 1
    else
        log_message "Successfully created: $file ($(du -h $file | cut -f1))"
    fi
}

# Helper function to verify required commands are available
check_command() {
    if ! command -v $1 &> /dev/null; then
        log_message "ERROR: $1 command not found. Please ensure it is installed and in PATH"
        exit 1
    fi
}

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Verify MultiQC is available
log_message "Checking required commands..."
check_command multiqc
check_command jq

# Define the output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/4_multiqc"

# Set and change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create output directories for MultiQC reports and QC data
log_message "Creating output directories..."
mkdir -p ${OUTPUT_DIR} || { log_message "ERROR: Failed to create output directory"; exit 1; }

# Run MultiQC to aggregate QC metrics
log_message "Running MultiQC to aggregate QC metrics..."
multiqc \
    --filename "multiqc_final" \
    --outdir ${OUTPUT_DIR} \
    --title "Final Sample QC Report" \
    --module fastqc \
    --module picard \
    --module samtools \
    --ignore "*_R1_001_fastqc*" \
    --ignore "*_R2_001_fastqc*" \
    --ignore "*_R1_paired_fastqc*" \
    --ignore "*_R2_paired_fastqc*" \
    --ignore "*_prefilter_flagstat*" \
    --ignore "*_prededup_flagstat*" \
    analysis/3_alignment/qc/*_flagstat.txt \
    analysis/3_alignment/qc/*_dup_metrics.txt \
    analysis/3_alignment/qc/*_insert_size_metrics.txt \
    2> logs/multiqc_final.log

# Create QC summary CSV with header
echo "Sample,RawReads,CleanReads,TrimRetentionRate,MappedReads,UniqueReads,FinalRetentionRate" > ${OUTPUT_DIR}/qc_summary.csv

# Process each sample to calculate QC metrics
for sample in GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3; do
    log_message "Processing statistics for ${sample}..."
    
    # Read trimming metrics from JSON
    if [[ -f analysis/2_quality_control/qc/${sample}_trimming_metrics.json ]]; then
        trim_metrics=$(cat analysis/2_quality_control/qc/${sample}_trimming_metrics.json)
        raw_reads=$(echo $trim_metrics | jq -r '.raw_r1_reads')
        clean_reads=$(echo $trim_metrics | jq -r '.clean_r1_reads')
        trim_retention=$(echo $trim_metrics | jq -r '.retention_rate')
    else
        log_message "WARNING: Trimming metrics not found for ${sample}"
        continue
    fi
    
    # Read alignment metrics from summary
    if [[ -f analysis/3_alignment/qc/${sample}_alignment_summary.txt ]]; then
        mapped_reads=$(grep "Mapped reads:" analysis/3_alignment/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' ')
        unique_reads=$(grep "Final reads after duplicate marking:" analysis/3_alignment/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' ')
        final_retention=$(grep "Overall retention rate:" analysis/3_alignment/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' %')
    else
        log_message "WARNING: Alignment metrics not found for ${sample}"
        continue
    fi
    
    # Log statistics for this sample
    log_message "${sample} statistics:"
    log_message "  Raw reads: ${raw_reads}"
    log_message "  Clean reads: ${clean_reads}"
    log_message "  Trim retention rate: ${trim_retention}%"
    log_message "  Mapped reads: ${mapped_reads}"
    log_message "  Unique reads: ${unique_reads}"
    log_message "  Final retention rate: ${final_retention}%"
    
    # Add statistics to summary CSV
    echo "${sample},${raw_reads},${clean_reads},${trim_retention},${mapped_reads},${unique_reads},${final_retention}" >> ${OUTPUT_DIR}/qc_summary.csv
done

check_output "${OUTPUT_DIR}/qc_summary.csv"
check_output "${OUTPUT_DIR}/multiqc_final.html"

log_message "MultiQC report generation completed successfully"