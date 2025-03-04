#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --output=logs/download_sra_%A_%a.out
#SBATCH --error=logs/download_sra_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=0-5

# Script to download SRA files for SRF_SES_V5 project
# These files are from a study on glioblastoma cells (SNB19) with SES and SOX2 transduction

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/data_from_ncbi_corrected

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create a log directory
mkdir -p logs

# Check if SRA toolkit is installed
if ! command -v prefetch &> /dev/null; then
    echo "SRA toolkit is not installed or not in PATH. Please install it first."
    echo "Visit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit"
    exit 1
fi

# Array of SRR accessions to download
declare -a SRR_ACCESSIONS=(
    "SRR18590288"  # SES-transduced SNB19
    "SRR18590287"  # SES-transduced SNB19
    "SRR18590286"  # SES-transduced SNB19
    "SRR18590285"  # SOX2-transduced SNB19
    "SRR18590284"  # SOX2-transduced SNB19
    "SRR18590283"  # SOX2-transduced SNB19
)

# Create a metadata file with sample information (only in the first job)
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    cat > sample_metadata.tsv << EOF
SRR_Accession	BioSample	GSM_ID	Cell_Line	Treatment
SRR18590288	SAMN27279093	GSM6008245	SNB19	SES-transduced
SRR18590287	SAMN27279092	GSM6008246	SNB19	SES-transduced
SRR18590286	SAMN27279091	GSM6008247	SNB19	SES-transduced
SRR18590285	SAMN27279090	GSM6008248	SNB19	SOX2-transduced
SRR18590284	SAMN27279089	GSM6008249	SNB19	SOX2-transduced
SRR18590283	SAMN27279088	GSM6008250	SNB19	SOX2-transduced
EOF
fi

# Get the current SRR accession based on the array task ID
current_accession=${SRR_ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

echo "Starting download of SRA file: $current_accession"
echo "This may take a while depending on your internet connection."
echo "Files will be downloaded to the current directory."

# Download the current SRR file
echo "Downloading $current_accession..."

# Use prefetch to download the SRA file
prefetch "$current_accession" --output-directory ./

if [ $? -eq 0 ]; then
    echo "Successfully downloaded $current_accession"
    
    # Convert SRA to fastq (optional - uncomment if needed)
    # echo "Converting $current_accession to FASTQ format..."
    # fasterq-dump "$current_accession" --split-files --outdir ./
else
    echo "Error downloading $current_accession."
fi

echo "Download process completed for $current_accession." 