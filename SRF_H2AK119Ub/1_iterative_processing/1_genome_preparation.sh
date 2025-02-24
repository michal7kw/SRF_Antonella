#!/bin/bash
#SBATCH --job-name=1_genome_preparation
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/1_genome_preparation.err"
#SBATCH --output="logs/1_genome_preparation.out"

# Description:
# This script automates the preparation of the human reference genome (GRCh38)
# for use in downstream ChIP-seq analysis pipelines. It downloads the genome
# sequence, builds a Bowtie2 index for read alignment, and generates a chromosome
# sizes file required by various genome analysis tools.

# Workflow:
# 1. Create the necessary directory structure to store the genome and log files.
# 2. Download the primary assembly of the human genome (GRCh38) from Ensembl.
#    The primary assembly excludes alternative contigs, ensuring a consistent reference.
# 3. Build Bowtie2 indices from the downloaded genome FASTA file. These indices
#    are essential for efficient alignment of sequencing reads to the genome.
# 4. Generate a chromosome sizes file using samtools and cut. This file lists
#    each chromosome and its length, which is used by tools for normalization
#    and visualization.

# Input Files:
#   - None. The script downloads the reference genome directly from the Ensembl FTP server.

# Output Files:
#   - $COMMON_DATA_DIR/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa:
#     The downloaded human reference genome FASTA file (GRCh38 primary assembly).
#   - $COMMON_DATA_DIR/genome/GRCh38.*.bt2:
#     Bowtie2 index files. These files are used by Bowtie2 to efficiently align
#     sequencing reads to the reference genome. Several files with the .bt2 extension
#     are created.
#   - $COMMON_DATA_DIR/genome/hg38.chrom.sizes:
#     A tab-separated file containing chromosome names and their corresponding lengths.
#     This file is used by various downstream analysis tools.
#   - logs/1_genome_preparation.out:
#     Standard output log file, capturing the script's progress and any messages.
#   - logs/1_genome_preparation.err:
#     Error log file, recording any errors or warnings encountered during script execution.

# Directory Structure:
#   - $COMMON_DATA_DIR/genome/:
#     This directory stores the downloaded reference genome FASTA file and the
#     generated Bowtie2 index files.
#   - logs/:
#     This directory stores the standard output and error log files for the script.

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

COMMON_DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA"

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

log_message "Creating necessary directories..."
# Create necessary directories
mkdir -p $COMMON_DATA_DIR/genome 
mkdir -p logs

# Download and prepare reference genome
cd $COMMON_DATA_DIR/genome

log_message "Downloading reference genome..."
# Download primary assembly to avoid alternative contigs
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build bowtie2 indices
log_message "Building human genome index..."
bowtie2-build --threads 32 Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38

# Generate chromosome sizes file for downstream analysis
log_message "Generating chromosome sizes file..."
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > hg38.chrom.sizes

log_message "Genome preparation completed"