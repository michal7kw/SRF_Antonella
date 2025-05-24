#!/bin/bash

# Restart processing for a specific failed sample
# Usage: ./restart_failed_sample.sh GSM6008246_SRR18590287

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <sample_name>"
    echo "Example: $0 GSM6008246_SRR18590287"
    echo ""
    echo "Available samples:"
    ls raw_data/*_R1.fastq.gz | sed 's/raw_data\///g' | sed 's/_R1.fastq.gz//g'
    exit 1
fi

SAMPLE_NAME=$1
R1_FILE="raw_data/${SAMPLE_NAME}_R1.fastq.gz"

if [[ ! -f "$R1_FILE" ]]; then
    echo "ERROR: Sample file not found: $R1_FILE"
    exit 1
fi

echo "🔄 Restarting processing for: $SAMPLE_NAME"
echo "=================================="

# Set up environment variables
export USE_SAMTOOLS_MARKDUP=true
export SKIP_DEDUP=false
export MAX_JOBS=1
export THREADS_PER_JOB=8  # Use more threads for single sample

# Source the process_sample function from the main script
source <(grep -A 200 "^process_sample()" create_big_wig_parallel_fixed.sh | head -n -1)

# Clean up any partial files
echo "Cleaning up partial files..."
rm -f aligned_bam/${SAMPLE_NAME}.*
rm -f bigwig_files/${SAMPLE_NAME}.*
rm -f logs/${SAMPLE_NAME}_*
rm -f qc_reports/${SAMPLE_NAME}_*

# Process the sample
echo "Starting processing..."
if process_sample "$R1_FILE"; then
    echo "✅ Successfully reprocessed $SAMPLE_NAME"
    echo ""
    echo "Output files:"
    echo "  BAM: aligned_bam/${SAMPLE_NAME}.dedup.bam"
    echo "  BigWig: bigwig_files/${SAMPLE_NAME}.bw"
    echo "  BigWig CPM: bigwig_files/${SAMPLE_NAME}.CPM.bw"
else
    echo "❌ Failed to reprocess $SAMPLE_NAME"
    echo "Check logs in logs/${SAMPLE_NAME}_*.log"
fi 