#!/bin/bash

# Quick-start wrapper for the parallelized Cut&Tag pipeline
# Automatically handles Java compatibility and runs the pipeline

echo "🚀 Cut&Tag Pipeline Quick Start"
echo "==============================="

# Check if we're in the right directory
if [[ ! -d "raw_data" ]]; then
    echo "❌ ERROR: raw_data directory not found"
    echo "   Please run this script from the SRF_SES_V5 directory"
    exit 1
fi

# Count fastq files
R1_COUNT=$(ls raw_data/*_R1.fastq.gz 2>/dev/null | wc -l)
echo "📁 Found $R1_COUNT samples in raw_data/"

if [[ $R1_COUNT -eq 0 ]]; then
    echo "❌ No R1 fastq.gz files found in raw_data/"
    exit 1
fi

# Check Java version and set default for non-interactive mode
if command -v java &> /dev/null; then
    JAVA_VERSION=$(java -version 2>&1 | grep -oP 'version "([0-9]+)' | grep -oP '[0-9]+' | head -1)
    echo "☕ Java version: $JAVA_VERSION"
    
    if [[ $JAVA_VERSION -lt 17 ]]; then
        echo "⚠️  Java $JAVA_VERSION detected (Picard needs 17+)"
        echo "✅ Will use samtools markdup instead (recommended)"
        export USE_SAMTOOLS_MARKDUP=true
        export SKIP_DEDUP=false
    else
        echo "✅ Java $JAVA_VERSION is compatible with Picard"
        export USE_SAMTOOLS_MARKDUP=false
        export SKIP_DEDUP=false
    fi
else
    echo "❌ Java not found! Please install Java:"
    echo "   conda install openjdk=17"
    exit 1
fi

# Set reasonable defaults for 20-core system
export MAX_JOBS=5
export THREADS_PER_JOB=4

echo ""
echo "🔧 Configuration:"
echo "   Parallel jobs: $MAX_JOBS"
echo "   Threads per job: $THREADS_PER_JOB" 
echo "   Total cores: $((MAX_JOBS * THREADS_PER_JOB))"
echo "   Duplicate removal: $([ "$USE_SAMTOOLS_MARKDUP" = "true" ] && echo "samtools markdup" || echo "Picard")"
echo ""

# Ask for confirmation
read -p "🚀 Start processing? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Cancelled by user"
    exit 1
fi

# Run the pipeline
echo "🏃 Starting pipeline..."
echo "==============================="

if [[ -f "create_big_wig_parallel_fixed.sh" ]]; then
    chmod +x create_big_wig_parallel_fixed.sh
    ./create_big_wig_parallel_fixed.sh
else
    echo "❌ Fixed pipeline script not found!"
    echo "   Expected: create_big_wig_parallel_fixed.sh"
    exit 1
fi 