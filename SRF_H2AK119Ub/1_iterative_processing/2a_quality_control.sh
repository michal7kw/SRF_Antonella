#!/bin/bash
#SBATCH --job-name=2_quality_control
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_quality_control_%a.err"
#SBATCH --output="logs/2_quality_control_%a.out"
#SBATCH --array=0-5  # Adjusted for all samples (6 samples total)

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
check_command fastqc
check_command trimmomatic

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Define samples
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/{fastqc/{pre_trim,post_trim},trimmed,qc} logs/trimming || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check if input files exist
if [[ ! -f ../DATA/fastq/${sample}_R1_001.fastq.gz ]] || [[ ! -f ../DATA/fastq/${sample}_R2_001.fastq.gz ]]; then
    log_message "ERROR: Input fastq files not found for ${sample}"
    log_message "Expected files:"
    log_message "  - ../DATA/fastq/${sample}_R1_001.fastq.gz"
    log_message "  - ../DATA/fastq/${sample}_R2_001.fastq.gz"
    exit 1
fi

log_message "Processing sample: ${sample}"

# Initial FastQC
log_message "Running initial FastQC for ${sample}..."
fastqc -o analysis/fastqc/pre_trim -t 16 \
    ../DATA/fastq/${sample}_R1_001.fastq.gz \
    ../DATA/fastq/${sample}_R2_001.fastq.gz

# Check FastQC output
check_output "analysis/fastqc/pre_trim/${sample}_R1_001_fastqc.html"
check_output "analysis/fastqc/pre_trim/${sample}_R2_001_fastqc.html"

# CUT&Tag-specific trimming parameters
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

# Check trimming output
check_output "analysis/trimmed/${sample}_R1_paired.fastq.gz"
check_output "analysis/trimmed/${sample}_R2_paired.fastq.gz"
check_output "logs/trimming/${sample}.log"

# Post-trimming FastQC
log_message "Running post-trim FastQC for ${sample}..."
fastqc -o analysis/fastqc/post_trim -t 16 \
    analysis/trimmed/${sample}_R1_paired.fastq.gz \
    analysis/trimmed/${sample}_R2_paired.fastq.gz

# Check post-trim FastQC output
check_output "analysis/fastqc/post_trim/${sample}_R1_paired_fastqc.html"
check_output "analysis/fastqc/post_trim/${sample}_R2_paired_fastqc.html"

# Calculate and log basic statistics
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

log_message "Quality control completed for ${sample}"