#!/bin/bash
#SBATCH --job-name=sra_to_fastq
#SBATCH --output=logs/sra_to_fastq_%A_%a.out
#SBATCH --error=logs/sra_to_fastq_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=0-5

# Script to convert SRA files to FASTQ format for SRF_SES_V5 project
# These files are from a study on glioblastoma cells (SNB19) with SES and SOX2 transduction

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/data_from_ncbi_corrected

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create a log directory
mkdir -p logs
# Create a directory for FASTQ files
mkdir -p fastq

# Check if SRA toolkit is installed
if ! command -v fasterq-dump &> /dev/null; then
    echo "SRA toolkit is not installed or not in PATH. Please install it first."
    echo "Visit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit"
    exit 1
fi

# Array of SRR accessions to process
declare -a SRR_ACCESSIONS=(
    "SRR18590288"  # SES-transduced SNB19
    "SRR18590287"  # SES-transduced SNB19
    "SRR18590286"  # SES-transduced SNB19
    "SRR18590285"  # SOX2-transduced SNB19
    "SRR18590284"  # SOX2-transduced SNB19
    "SRR18590283"  # SOX2-transduced SNB19
)

# Get the current SRR accession based on the array task ID
current_accession=${SRR_ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

echo "Starting conversion of SRA file to FASTQ: $current_accession"

# Check if the SRA file exists
if [ ! -f "./$current_accession/$current_accession.sra" ]; then
    echo "SRA file for $current_accession not found. Please download it first."
    exit 1
fi

# Convert SRA to FASTQ
echo "Converting $current_accession to FASTQ format..."
fasterq-dump "$current_accession/$current_accession.sra" --split-files --threads 4 --outdir ./fastq

if [ $? -eq 0 ]; then
    echo "Successfully converted $current_accession to FASTQ format"
    
    # Compress the FASTQ files to save space
    echo "Compressing FASTQ files..."
    gzip ./fastq/${current_accession}_*.fastq
    
    # Create symbolic links with more descriptive names based on metadata
    treatment=$(grep "$current_accession" sample_metadata.tsv | cut -f5)
    replicate=$(grep -n "$current_accession" sample_metadata.tsv | grep -o "SES-transduced" | wc -l)
    
    if [[ "$treatment" == "SES-transduced" ]]; then
        # Find which replicate this is among SES samples
        ses_samples=($(grep "SES-transduced" sample_metadata.tsv | cut -f1))
        for i in "${!ses_samples[@]}"; do
            if [[ "${ses_samples[$i]}" == "$current_accession" ]]; then
                replicate=$((i+1))
                break
            fi
        done
        ln -sf ./fastq/${current_accession}_1.fastq.gz ./fastq/SES_rep${replicate}_R1.fastq.gz
        ln -sf ./fastq/${current_accession}_2.fastq.gz ./fastq/SES_rep${replicate}_R2.fastq.gz
    elif [[ "$treatment" == "SOX2-transduced" ]]; then
        # Find which replicate this is among SOX2 samples
        sox2_samples=($(grep "SOX2-transduced" sample_metadata.tsv | cut -f1))
        for i in "${!sox2_samples[@]}"; do
            if [[ "${sox2_samples[$i]}" == "$current_accession" ]]; then
                replicate=$((i+1))
                break
            fi
        done
        ln -sf ./fastq/${current_accession}_1.fastq.gz ./fastq/SOX2_rep${replicate}_R1.fastq.gz
        ln -sf ./fastq/${current_accession}_2.fastq.gz ./fastq/SOX2_rep${replicate}_R2.fastq.gz
    fi
    
    echo "Created symbolic links with descriptive names for $current_accession"
else
    echo "Error converting $current_accession to FASTQ format."
fi

echo "Conversion process completed for $current_accession." 