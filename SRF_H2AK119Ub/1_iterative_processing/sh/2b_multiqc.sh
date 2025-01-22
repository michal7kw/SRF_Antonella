#!/bin/bash
#SBATCH --job-name=2_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_multiqc.err"
#SBATCH --output="logs/2_multiqc.out"

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
check_command multiqc

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/{multiqc,qc} || { log_message "ERROR: Failed to create directories"; exit 1; }

# Check if FastQC results exist
if [[ ! -d analysis/fastqc/pre_trim ]] || [[ ! -d analysis/fastqc/post_trim ]]; then
    log_message "ERROR: FastQC directories not found. Please run quality control first."
    exit 1
fi

# Run MultiQC for pre-trimming data
log_message "Creating pre-trimming MultiQC report..."
multiqc \
    --filename "multiqc_pre_trimming" \
    --outdir analysis/multiqc \
    --title "Pre-trimming QC Report" \
    analysis/fastqc/pre_trim/*_fastqc* \
    2> logs/multiqc_pre_trim.log

check_output "analysis/multiqc/multiqc_pre_trimming.html"

# Run MultiQC for post-trimming data
log_message "Creating post-trimming MultiQC report..."
multiqc \
    --filename "multiqc_post_trimming" \
    --outdir analysis/multiqc \
    --title "Post-trimming QC Report" \
    analysis/fastqc/post_trim/*_paired_fastqc* \
    2> logs/multiqc_post_trim.log

check_output "analysis/multiqc/multiqc_post_trimming.html"

# Generate QC summary
log_message "Generating QC summary..."
echo "Sample,RawReads,CleanReads,RetentionRate" > analysis/qc/qc_summary.csv

for sample in GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3; do
    log_message "Processing statistics for ${sample}..."
    
    # Check if input files exist
    if [[ ! -f ../DATA/fastq/${sample}_R1_001.fastq.gz ]]; then
        log_message "WARNING: Raw fastq file not found for ${sample}, skipping..."
        continue
    fi
    
    if [[ ! -f analysis/trimmed/${sample}_R1_paired.fastq.gz ]]; then
        log_message "WARNING: Cleaned fastq file not found for ${sample}, skipping..."
        continue
    fi
    
    # Calculate statistics
    raw_reads=$(zcat ../DATA/fastq/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))
    clean_reads=$(zcat analysis/trimmed/${sample}_R1_paired.fastq.gz | echo $((`wc -l`/4)))
    retention_rate=$(echo "scale=2; $clean_reads/$raw_reads * 100" | bc)
    
    log_message "${sample} statistics:"
    log_message "  Raw reads: ${raw_reads}"
    log_message "  Clean reads: ${clean_reads}"
    log_message "  Retention rate: ${retention_rate}%"
    
    echo "${sample},${raw_reads},${clean_reads},${retention_rate}" >> analysis/qc/qc_summary.csv
done

check_output "analysis/qc/qc_summary.csv"

log_message "MultiQC report generation completed successfully"