#!/bin/bash

# This script generates a chromosome sizes file from a genome FASTA file

GENOME_FASTA="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTPUT_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes"

echo "Generating chromosome sizes file from FASTA..."
samtools faidx $GENOME_FASTA
cut -f1,2 ${GENOME_FASTA}.fai > $OUTPUT_FILE
echo "Chromosome sizes file created at: $OUTPUT_FILE"
