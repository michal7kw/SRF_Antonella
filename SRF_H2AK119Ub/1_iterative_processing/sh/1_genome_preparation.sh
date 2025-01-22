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

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

log_message "Creating necessary directories..."
# Create necessary directories
mkdir -p genome logs
mkdir -p analysis/{fastqc/{pre_trim,post_trim},trimmed,aligned,peaks,visualization,qc,annotation}

# Download and prepare reference genome
cd genome

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

cd ..

log_message "Genome preparation completed"