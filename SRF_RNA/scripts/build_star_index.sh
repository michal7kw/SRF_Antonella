#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --output=logs/build_index.out
#SBATCH --error=logs/build_index.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq

# Exit on error
set -e

# Load conda environment (if needed)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Parameters
GENOME_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome_star"
GENOME_FASTA="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.nochr.gtf"
THREADS=32

# Create output directory and clean it
echo "Creating and cleaning genome directory..."
rm -rf ${GENOME_DIR}/*
mkdir -p ${GENOME_DIR}

# Build STAR index
echo "Building STAR index..."
STAR --runMode genomeGenerate \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${GENOME_FASTA} \
     --sjdbGTFfile ${GTF_FILE} \
     --runThreadN ${THREADS} \
     --sjdbOverhang 100 \
     --limitGenomeGenerateRAM 60000000000

# Create flag file
echo "Creating flag file..."
touch /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/logs/star_index/index_complete.flag

echo "STAR index generation complete."
