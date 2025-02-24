#!/bin/bash
#SBATCH --job-name=2_quality_control
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_quality_control.err"
#SBATCH --output="logs/2_quality_control.out"
#SBATCH --array=0-5  # Process 6 samples in parallel (GFP_1-3 and YAF_1-3)

### Optional ###
# SBATCH --exclusive

# Hardcode the step to metrics
# STEP="metrics"

# Get the step parameter first, before any other operations
STEP="${1:-all}"

# Documentation:
# Quality Control and Trimming Pipeline for CUT&Tag Data
# This script performs comprehensive quality control and read processing for paired-end sequencing data
# from CUT&Tag experiments. It processes multiple samples in parallel using a SLURM array job to ensure
# efficient and reliable data processing. The pipeline includes initial quality assessment, adapter trimming,
# quality filtering, and post-trimming quality assessment. It generates detailed read statistics and QC metrics
# in JSON format for downstream analysis and reporting.

# Key Features:
# - Parallel processing of 6 samples (GFP_1-3 and YAF_1-3) using SLURM array jobs.
# - Initial quality assessment of raw reads using FastQC.
# - Adapter trimming and quality filtering using Trimmomatic with specific parameters for CUT&Tag data.
# - Post-trimming quality assessment using FastQC to evaluate the effect of trimming.
# - Calculation of read statistics (raw reads, clean reads, retention rate) to quantify the trimming process.
# - Generation of comprehensive QC reports, including FastQC reports, Trimmomatic logs, and JSON-formatted metrics.

# Pipeline Steps:
# 1. Initial quality assessment with FastQC on raw reads to identify potential issues.
# 2. Read trimming and quality filtering with Trimmomatic to remove adapters and low-quality bases.
# 3. Post-trimming quality assessment with FastQC to verify the quality of trimmed reads.
# 4. Calculation of read statistics and retention rates to quantify the number of reads retained after trimming.
# 5. Generation of comprehensive QC reports in HTML (FastQC) and JSON formats for easy interpretation and integration with downstream analysis.

# Required Tools:
# - FastQC (v0.11.9 or later): For assessing the quality of raw and trimmed reads.
# - Trimmomatic (v0.39 or later): For trimming adapters and filtering low-quality reads.
# - GNU core utilities (wc, zcat, etc.): For basic file manipulation and counting reads.
# - bc (for calculations): For calculating retention rates.
# - conda: For environment and package management

# Input Files (per sample):
# - ../DATA/fastq/{sample}_R1_001.fastq.gz: Raw forward reads in FASTQ format (gzip compressed).
# - ../DATA/fastq/{sample}_R2_001.fastq.gz: Raw reverse reads in FASTQ format (gzip compressed).
# - ../DATA/TruSeq3-PE.fa: FASTA file containing adapter sequences for Trimmomatic.

# Output Files (per sample):
# Quality Control Reports:
# - ${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/{sample}_R1_001_fastqc.html: Initial FastQC report for forward reads (R1).
# - ${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/{sample}_R1_001_fastqc.zip: Initial FastQC zip archive for forward reads (R1).
# - ${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/{sample}_R2_001_fastqc.html: Initial FastQC report for reverse reads (R2).
# - ${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/{sample}_R2_001_fastqc.zip: Initial FastQC zip archive for reverse reads (R2).
# - ${OUTPUT_DIR}/2_trimming/fastqc/post_trim/{sample}_R1_paired_fastqc.html: Post-trim FastQC report for paired forward reads (R1).
# - ${OUTPUT_DIR}/2_trimming/fastqc/post_trim/{sample}_R1_paired_fastqc.zip: Post-trim FastQC zip archive for paired forward reads (R1).
# - ${OUTPUT_DIR}/2_trimming/fastqc/post_trim/{sample}_R2_paired_fastqc.html: Post-trim FastQC report for paired reverse reads (R2).
# - ${OUTPUT_DIR}/2_trimming/fastqc/post_trim/{sample}_R2_paired_fastqc.zip: Post-trim FastQC zip archive for paired reverse reads (R2).

# Trimmed Reads:
# - ${OUTPUT_DIR}/2_trimming/{sample}_R1_paired.fastq.gz: Trimmed and paired forward reads in FASTQ format (gzip compressed).
# - ${OUTPUT_DIR}/2_trimming/{sample}_R1_unpaired.fastq.gz: Unpaired forward reads in FASTQ format (gzip compressed).
# - ${OUTPUT_DIR}/2_trimming/{sample}_R2_paired.fastq.gz: Trimmed and paired reverse reads in FASTQ format (gzip compressed).
# - ${OUTPUT_DIR}/2_trimming/{sample}_R2_unpaired.fastq.gz: Unpaired reverse reads in FASTQ format (gzip compressed).

# Logs and Metrics:
# - ${OUTPUT_DIR}/2_trimming/logs/trimming/{sample}.log: Detailed Trimmomatic log file containing trimming statistics.
# - ${OUTPUT_DIR}/2_trimming/qc/{sample}_trimming_metrics.json: JSON file containing read statistics and metrics (raw reads, clean reads, retention rate).
# - logs/2_quality_control.out: Standard output log for the entire script.
# - logs/2_quality_control.err: Error log for the entire script.

# Exit script on error, undefined variables, and pipe failures
set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check if output exists and is not empty
# Used to verify successful generation of output files
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
# Used to verify required tools are available
check_command() {
    if ! command -v $1 &> /dev/null; then
        log_message "ERROR: $1 command not found. Please ensure it is installed and in PATH"
        exit 1
    fi
}

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Activate the specific environment before processing step parameter

# Verify required tools are available
log_message "Checking required commands..."
check_command fastqc
check_command trimmomatic

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis"

# Define sample names array and get current sample based on array task ID
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Create directory structure for output files
log_message "Creating output directories..."
mkdir -p ${OUTPUT_DIR}/{2_trimming/{fastqc/{pre_trim,post_trim},logs/trimming,qc},2_trimming} logs/trimming || { log_message "ERROR: Failed to create directories"; exit 1; }

# Verify input files exist
if [[ ! -f ../DATA/fastq/${sample}_R1_001.fastq.gz ]] || [[ ! -f ../DATA/fastq/${sample}_R2_001.fastq.gz ]]; then
    log_message "ERROR: Input fastq files not found for ${sample}"
    log_message "Expected files:"
    log_message "  - ../DATA/fastq/${sample}_R1_001.fastq.gz"
    log_message "  - ../DATA/fastq/${sample}_R2_001.fastq.gz"
    exit 1
fi

log_message "Processing sample: ${sample}"

# Function to run initial FastQC
run_initial_fastqc() {
    local sample=$1
    log_message "Running initial FastQC for ${sample}..."
    fastqc -o ${OUTPUT_DIR}/2_trimming/fastqc/pre_trim -t 16 \
        ../DATA/fastq/${sample}_R1_001.fastq.gz \
        ../DATA/fastq/${sample}_R2_001.fastq.gz

    check_output "${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/${sample}_R1_001_fastqc.html"
    check_output "${OUTPUT_DIR}/2_trimming/fastqc/pre_trim/${sample}_R2_001_fastqc.html"
}

# Function to run Trimmomatic
run_trimmomatic() {
    local sample=$1
    log_message "Trimming ${sample}..."
    trimmomatic PE -threads 16 \
        ../DATA/fastq/${sample}_R1_001.fastq.gz \
        ../DATA/fastq/${sample}_R2_001.fastq.gz \
        ${OUTPUT_DIR}/2_trimming/${sample}_R1_paired.fastq.gz \
        ${OUTPUT_DIR}/2_trimming/${sample}_R1_unpaired.fastq.gz \
        ${OUTPUT_DIR}/2_trimming/${sample}_R2_paired.fastq.gz \
        ${OUTPUT_DIR}/2_trimming/${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:../DATA/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 TRAILING:20 \
        SLIDINGWINDOW:4:15 \
        MINLEN:25 \
        2> ${OUTPUT_DIR}/2_trimming/logs/trimming/${sample}.log

    check_output "${OUTPUT_DIR}/2_trimming/${sample}_R1_paired.fastq.gz"
    check_output "${OUTPUT_DIR}/2_trimming/${sample}_R2_paired.fastq.gz"
    check_output "${OUTPUT_DIR}/2_trimming/logs/trimming/${sample}.log"
}

# Function to run post-trim FastQC
run_post_trim_fastqc() {
    local sample=$1
    log_message "Running post-trim FastQC for ${sample}..."
    fastqc -o ${OUTPUT_DIR}/2_trimming/fastqc/post_trim -t 16 \
        ${OUTPUT_DIR}/2_trimming/${sample}_R1_paired.fastq.gz \
        ${OUTPUT_DIR}/2_trimming/${sample}_R2_paired.fastq.gz

    check_output "${OUTPUT_DIR}/2_trimming/fastqc/post_trim/${sample}_R1_paired_fastqc.html"
    check_output "${OUTPUT_DIR}/2_trimming/fastqc/post_trim/${sample}_R2_paired_fastqc.html"
}

# Function to generate QC metrics
generate_qc_metrics() {
    local sample=$1
    log_message "Calculating statistics for ${sample}..."
    raw_r1_reads=$(zcat ../DATA/fastq/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))
    raw_r2_reads=$(zcat ../DATA/fastq/${sample}_R2_001.fastq.gz | echo $((`wc -l`/4)))
    clean_r1_reads=$(zcat ${OUTPUT_DIR}/2_trimming/${sample}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
    clean_r2_reads=$(zcat ${OUTPUT_DIR}/2_trimming/${sample}_R2_paired.fastq.gz | echo $((`wc -l`/4)))

    log_message "Statistics for ${sample}:"
    log_message "  Raw R1 reads: ${raw_r1_reads}"
    log_message "  Raw R2 reads: ${raw_r2_reads}"
    log_message "  Clean R1 reads: ${clean_r1_reads}"
    log_message "  Clean R2 reads: ${clean_r2_reads}"
    log_message "  Retention rate: $(echo "scale=2; $clean_r1_reads/$raw_r1_reads * 100" | bc)%"

    # Store detailed QC metrics in JSON format
    cat << EOF > ${OUTPUT_DIR}/2_trimming/qc/${sample}_trimming_metrics.json
{
    "sample": "${sample}",
    "raw_r1_reads": ${raw_r1_reads},
    "raw_r2_reads": ${raw_r2_reads},
    "clean_r1_reads": ${clean_r1_reads},
    "clean_r2_reads": ${clean_r2_reads},
    "retention_rate": $(echo "scale=2; $clean_r1_reads/$raw_r1_reads * 100" | bc)
}
EOF
}

# Main execution based on step parameter
case $STEP in
    "fastqc")
        run_initial_fastqc ${sample}
        ;;
    "trim")
        run_trimmomatic ${sample}
        ;;
    "post_fastqc")
        run_post_trim_fastqc ${sample}
        ;;
    "metrics")
        generate_qc_metrics ${sample}
        ;;
    "all")
        run_initial_fastqc ${sample}
        run_trimmomatic ${sample}
        run_post_trim_fastqc ${sample}
        generate_qc_metrics ${sample}
        ;;
    *)
        log_message "ERROR: Invalid step specified. Valid steps are: fastqc, trim, post_fastqc, metrics, all"
        exit 1
        ;;
esac

log_message "Quality control completed for ${sample}"