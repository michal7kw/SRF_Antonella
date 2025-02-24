#!/bin/bash
#SBATCH --job-name=RNA
#SBATCH --output=logs/RNA.out
#SBATCH --error=logs/RNA.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=120:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Exit on error
set -e

# Create necessary directories
mkdir -p logs/cluster_logs
mkdir -p DATA
mkdir -p results/fastqc

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA

# Check if input files exist
for sample in C1 C2 C3 GFP1 GFP2 GFP3 YAF1 YAF2 YAF3; do
    if [ ! -f "DATA/${sample}/${sample}_L001_R1_001.fastq.gz" ] || [ ! -f "DATA/${sample}/${sample}_L001_R2_001.fastq.gz" ]; then
        echo "Error: Input files for sample ${sample} not found in DATA/${sample}/"
        echo "Please ensure the following files exist:"
        echo "  - DATA/${sample}/${sample}_L001_R1_001.fastq.gz"
        echo "  - DATA/${sample}/${sample}_L001_R2_001.fastq.gz"
        exit 1
    fi
done

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Unlock snakemake working directory if necessary
# snakemake --unlock

# Run snakemake with dry-run first to validate
echo "Performing dry-run to validate workflow..."
snakemake --snakefile Snakefile --dry-run --rerun-incomplete

if [ $? -eq 0 ]; then
    echo "Dry-run successful, starting actual run..."
    # Run snakemake
    snakemake \
        --snakefile Snakefile \
        --executor slurm \
        --rerun-incomplete \
        --jobs 100 \
        --default-resources \
            slurm_partition=workq \
            mem_mb=32000 \
            runtime=1440 \
            threads=8 \
            nodes=1 \
        --jobscript slurm-jobscript.sh \
        --latency-wait 60 \
        --rerun-incomplete \
        --keep-going
else
    echo "Dry-run failed, please check the error messages above"
    exit 1
fi