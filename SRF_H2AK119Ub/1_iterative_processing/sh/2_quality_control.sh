#!/bin/bash
#SBATCH --job-name=2_quality_control
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_quality_control_%a.err"
#SBATCH --output="logs/2_quality_control_%a.out"
#SBATCH --array=0-5  # Process 6 samples in parallel (GFP_1-3 and YAF_1-3)

### Optional ###
# SBATCH --exclusive

# Documentation:
# This script performs quality control and trimming on paired-end sequencing data
# It processes multiple samples in parallel using a SLURM array job
# For each sample, it:
# 1. Runs initial FastQC on raw reads
# 2. Trims reads using Trimmomatic with CUT&Tag-specific parameters
# 3. Runs FastQC again on trimmed reads
# 4. Calculates and logs read statistics

# Input files (per sample):
# - ../DATA/fastq/{sample}_R1_001.fastq.gz - Raw forward reads
# - ../DATA/fastq/{sample}_R2_001.fastq.gz - Raw reverse reads

# Output files (per sample):
# - analysis/fastqc/pre_trim/{sample}_R[1,2]_001_fastqc.html - Initial FastQC reports
# - analysis/trimmed/{sample}_R[1,2]_paired.fastq.gz - Trimmed paired reads
# - analysis/trimmed/{sample}_R[1,2]_unpaired.fastq.gz - Trimmed unpaired reads
# - analysis/fastqc/post_trim/{sample}_R[1,2]_paired_fastqc.html - Post-trim FastQC reports
# - logs/trimming/{sample}.log - Trimmomatic log file

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
conda activate snakemake

# Verify required tools are available
log_message "Checking required commands..."
check_command fastqc
check_command trimmomatic

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define sample names array and get current sample based on array task ID
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Create directory structure for output files
log_message "Creating output directories..."
mkdir -p analysis/{fastqc/{pre_trim,post_trim},trimmed,qc} logs/trimming || { log_message "ERROR: Failed to create directories"; exit 1; }

# Verify input files exist
if [[ ! -f ../DATA/fastq/${sample}_R1_001.fastq.gz ]] || [[ ! -f ../DATA/fastq/${sample}_R2_001.fastq.gz ]]; then
    log_message "ERROR: Input fastq files not found for ${sample}"
    log_message "Expected files:"
    log_message "  - ../DATA/fastq/${sample}_R1_001.fastq.gz"
    log_message "  - ../DATA/fastq/${sample}_R2_001.fastq.gz"
    exit 1
fi

log_message "Processing sample: ${sample}"

# Run initial FastQC on raw reads
log_message "Running initial FastQC for ${sample}..."
fastqc -o analysis/fastqc/pre_trim -t 16 \
    ../DATA/fastq/${sample}_R1_001.fastq.gz \
    ../DATA/fastq/${sample}_R2_001.fastq.gz

# Verify FastQC output
check_output "analysis/fastqc/pre_trim/${sample}_R1_001_fastqc.html"
check_output "analysis/fastqc/pre_trim/${sample}_R2_001_fastqc.html"

# Run Trimmomatic with CUT&Tag-specific parameters:
# - ILLUMINACLIP: Remove Illumina adapters
# - LEADING/TRAILING: Remove low quality bases from start/end
# - SLIDINGWINDOW: Remove low quality regions
# - MINLEN: Remove too short reads
log_message "Trimming ${sample}..."
trimmomatic PE -threads 16 \
    ../DATA/fastq/${sample}_R1_001.fastq.gz \
    ../DATA/fastq/${sample}_R2_001.fastq.gz \
    analysis/trimmed/${sample}_R1_paired.fastq.gz \
    analysis/trimmed/${sample}_R1_unpaired.fastq.gz \
    analysis/trimmed/${sample}_R2_paired.fastq.gz \
    analysis/trimmed/${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25 \
    2> logs/trimming/${sample}.log

# Verify trimming output
check_output "analysis/trimmed/${sample}_R1_paired.fastq.gz"
check_output "analysis/trimmed/${sample}_R2_paired.fastq.gz"
check_output "logs/trimming/${sample}.log"

# Run FastQC on trimmed reads
log_message "Running post-trim FastQC for ${sample}..."
fastqc -o analysis/fastqc/post_trim -t 16 \
    analysis/trimmed/${sample}_R1_paired.fastq.gz \
    analysis/trimmed/${sample}_R2_paired.fastq.gz

# Verify post-trim FastQC output
check_output "analysis/fastqc/post_trim/${sample}_R1_paired_fastqc.html"
check_output "analysis/fastqc/post_trim/${sample}_R2_paired_fastqc.html"

# Calculate and log read statistics
log_message "Calculating statistics for ${sample}..."
raw_r1_reads=$(zcat ../DATA/fastq/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))
raw_r2_reads=$(zcat ../DATA/fastq/${sample}_R2_001.fastq.gz | echo $((`wc -l`/4)))
clean_r1_reads=$(zcat analysis/trimmed/${sample}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
clean_r2_reads=$(zcat analysis/trimmed/${sample}_R2_paired.fastq.gz | echo $((`wc -l`/4)))

log_message "Statistics for ${sample}:"
log_message "  Raw R1 reads: ${raw_r1_reads}"
log_message "  Raw R2 reads: ${raw_r2_reads}"
log_message "  Clean R1 reads: ${clean_r1_reads}"
log_message "  Clean R2 reads: ${clean_r2_reads}"
log_message "  Retention rate: $(echo "scale=2; $clean_r1_reads/$raw_r1_reads * 100" | bc)%"

# Store detailed QC metrics in JSON format
cat << EOF > analysis/qc/${sample}_trimming_metrics.json
{
    "sample": "${sample}",
    "raw_r1_reads": ${raw_r1_reads},
    "raw_r2_reads": ${raw_r2_reads},
    "clean_r1_reads": ${clean_r1_reads},
    "clean_r2_reads": ${clean_r2_reads},
    "retention_rate": $(echo "scale=2; $clean_r1_reads/$raw_r1_reads * 100" | bc)
}
EOF

log_message "Quality control completed for ${sample}"