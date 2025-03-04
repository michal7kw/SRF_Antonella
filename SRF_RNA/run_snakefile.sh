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
mkdir -p logs/deseq2
mkdir -p logs/bigwig
mkdir -p DATA
mkdir -p results/fastqc
mkdir -p results/bigwig
mkdir -p results/deseq2
mkdir -p scripts

# Generate chromosome sizes file if it doesn't exist
if [ ! -f "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes" ]; then
    echo "Generating chromosome sizes file..."
    ./scripts/generate_chrom_sizes.sh
fi

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA

# Check if STAR index exists and is complete
STAR_INDEX="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome_star"
if [ ! -d "$STAR_INDEX" ] || [ ! -f "$STAR_INDEX/genomeParameters.txt" ]; then
    echo "ERROR: STAR index not found or incomplete. Running the build_star_index.sh script first."
    sbatch scripts/build_star_index.sh
    echo "STAR index build job submitted. Please wait for it to complete before running this script again."
    exit 1
fi

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

# Install required R packages if not already installed
# Rscript -e 'if(!require("DESeq2")) BiocManager::install("DESeq2", update=FALSE)'
# Rscript -e 'if(!require("tidyverse")) install.packages("tidyverse", repos="https://cloud.r-project.org")'
# Rscript -e 'if(!require("ggrepel")) install.packages("ggrepel", repos="https://cloud.r-project.org")'
# Rscript -e 'if(!require("pheatmap")) install.packages("pheatmap", repos="https://cloud.r-project.org")'
# Rscript -e 'if(!require("RColorBrewer")) install.packages("RColorBrewer", repos="https://cloud.r-project.org")'
# Rscript -e 'if(!require("EnhancedVolcano")) BiocManager::install("EnhancedVolcano", update=FALSE)'

# Unlock snakemake working directory if necessary
echo "Unlocking working directory..."
snakemake --unlock

# Run snakemake with dry-run first to validate
echo "Performing dry-run to validate workflow..."
# snakemake --snakefile Snakefile --dry-run --rerun-incomplete

if [ $? -eq 0 ]; then
    echo "Dry-run successful, starting actual run..."
    # Run snakemake
    snakemake \
        --snakefile Snakefile \
        --executor slurm \
        --jobs 100 \
        --default-resources \
            slurm_partition=workq \
            mem_mb=32000 \
            runtime=1440 \
            threads=8 \
            nodes=1 \
        --jobscript slurm-jobscript.sh \
        --latency-wait 60 \
        --keep-going
else
    echo "Dry-run failed, please check the error messages above"
    exit 1
fi

# --rerun-incomplete \
# --use-conda \
# --resources mem_mb=32000 \