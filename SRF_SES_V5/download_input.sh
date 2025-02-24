#!/bin/bash

# Create directory for the data if it doesn't exist
mkdir -p data_from_ncbi

# Download the SRA file
prefetch SRR18590303

# Convert SRA to fastq files
fasterq-dump --split-files SRR18590303 \
    --outdir data_from_ncbi \
    --threads 16

# Cleanup the SRA file to save space
rm -rf ~/ncbi/public/sra/SRR18590303.sra

# Update the INPUT_SAMPLE in Snakefile
sed -i 's/SRR18590290/SRR18590303/' Snakefile 