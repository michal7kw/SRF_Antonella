#!/bin/bash
#SBATCH --job-name=3b_multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3b_multiqc.err"
#SBATCH --output="logs/3b_multiqc.out"

# Documentation:
# This script runs MultiQC to aggregate FastQC reports and generate QC summaries
# It processes both pre- and post-trimming FastQC results and creates a QC metrics summary

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

# Set and change to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create output directories for MultiQC reports and QC data
log_message "Creating output directories..."
mkdir -p analysis/{multiqc,qc} || { log_message "ERROR: Failed to create directories"; exit 1; }

# # Verify FastQC results exist before proceeding
# if [[ ! -d analysis/fastqc/pre_trim ]] || [[ ! -d analysis/fastqc/post_trim ]]; then
#     log_message "ERROR: FastQC directories not found. Please run quality control first."
#     exit 1
# fi

# # Check if alignment results exist
# if [[ ! -f analysis/qc/*_alignment_summary.txt ]] || \
#    [[ ! -f analysis/qc/*_flagstat.txt ]] || \
#    [[ ! -f analysis/qc/*_dup_metrics.txt ]]; then
#     log_message "ERROR: Alignment results not found. Please run 3_alignment.sh first"
#     exit 1
# fi

# Modify the MultiQC command to include alignment metrics
multiqc \
    --filename "multiqc_final" \
    --outdir analysis/multiqc \
    --title "Complete Pipeline QC Report" \
    analysis/fastqc/pre_trim/*_fastqc* \
    analysis/fastqc/post_trim/*_paired_fastqc* \
    analysis/qc/*_flagstat.txt \
    analysis/qc/*_dup_metrics.txt \
    analysis/qc/*_insert_size_metrics.txt \
    2> logs/multiqc_final.log

# Modify the QC summary to include alignment metrics
echo "Sample,RawReads,CleanReads,TrimRetentionRate,MappedReads,UniqueReads,FinalRetentionRate" > analysis/qc/qc_summary.csv

# Process each sample to calculate QC metrics
for sample in GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3; do
    log_message "Processing statistics for ${sample}..."
    
    # # Verify input files exist
    # if [[ ! -f ../DATA/fastq/${sample}_R1_001.fastq.gz ]]; then
    #     log_message "WARNING: Raw fastq file not found for ${sample}, skipping..."
    #     continue
    # fi
    
    # if [[ ! -f analysis/trimmed/${sample}_R1_paired.fastq.gz ]]; then
    #     log_message "WARNING: Cleaned fastq file not found for ${sample}, skipping..."
    #     continue
    # fi
    
    # Read trimming metrics from JSON
    if [[ -f analysis/qc/${sample}_trimming_metrics.json ]]; then
        trim_metrics=$(cat analysis/qc/${sample}_trimming_metrics.json)
        raw_reads=$(echo $trim_metrics | jq -r '.raw_r1_reads')
        clean_reads=$(echo $trim_metrics | jq -r '.clean_r1_reads')
        trim_retention=$(echo $trim_metrics | jq -r '.retention_rate')
    else
        log_message "WARNING: Trimming metrics not found for ${sample}"
        continue
    fi
    
    # Read alignment metrics from summary
    if [[ -f analysis/qc/${sample}_alignment_summary.txt ]]; then
        mapped_reads=$(grep "Mapped reads:" analysis/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' ')
        unique_reads=$(grep "Final reads after duplicate marking:" analysis/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' ')
        final_retention=$(grep "Overall retention rate:" analysis/qc/${sample}_alignment_summary.txt | cut -d: -f2 | tr -d ' %')
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
    echo "${sample},${raw_reads},${clean_reads},${trim_retention},${mapped_reads},${unique_reads},${final_retention}" >> analysis/qc/qc_summary.csv
done

check_output "analysis/qc/qc_summary.csv"

log_message "MultiQC report generation completed successfully"