#!/bin/bash

# Cut&Tag to BigWig Processing Pipeline - Parallelized Version
# Processes all downloaded Cut&Tag data to create BigWig files using parallel processing

# Configuration
MAX_JOBS=5          # Number of samples to process in parallel (leave some cores for system)
THREADS_PER_JOB=4   # Threads per sample (5 jobs * 4 threads = 20 cores)
GENOME_INDEX="/home/michal/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as" 
CHROM_SIZES="/home/michal/COMMON_FILES/hg38.chrom.sizes"

# Create output directories
mkdir -p aligned_bam
mkdir -p bigwig_files
mkdir -p qc_reports
mkdir -p logs

# Function to process a single sample
process_sample() {
    local r1=$1
    local sample_base=$(basename $r1 _R1.fastq.gz)
    local gsm=$(echo $sample_base | cut -d'_' -f1)
    
    echo "========================================"
    echo "Processing: $sample_base (PID: $$)"
    echo "GSM: $gsm"
    echo "Threads: $THREADS_PER_JOB"
    echo "========================================"
    
    # Define file paths
    local r2=${r1/_R1.fastq.gz/_R2.fastq.gz}
    local bam="aligned_bam/${sample_base}.bam"
    local sorted_bam="aligned_bam/${sample_base}.sorted.bam"
    local dedup_bam="aligned_bam/${sample_base}.dedup.bam"
    local bigwig="bigwig_files/${sample_base}.bw"
    
    # Step 1: FastQC on raw reads
    echo "[$(date)] Running FastQC for $sample_base..."
    fastqc -t $THREADS_PER_JOB -o qc_reports $r1 $r2 2>&1 | tee logs/${sample_base}_fastqc.log
    
    # Step 2: Align with Bowtie2
    echo "[$(date)] Aligning $sample_base with Bowtie2..."
    bowtie2 -x $GENOME_INDEX \
        -1 $r1 \
        -2 $r2 \
        -p $THREADS_PER_JOB \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -I 10 -X 1000 \
        2> logs/${sample_base}_bowtie2.log \
        | samtools view -@ $THREADS_PER_JOB -bS - > $bam
    
    # Step 3: Sort BAM
    echo "[$(date)] Sorting BAM for $sample_base..."
    samtools sort -@ $THREADS_PER_JOB -o $sorted_bam $bam
    samtools index $sorted_bam
    
    # Remove unsorted BAM to save space
    rm $bam
    
    # Step 4: Mark and remove duplicates
    echo "[$(date)] Removing duplicates for $sample_base..."
    picard MarkDuplicates \
        I=$sorted_bam \
        O=$dedup_bam \
        M=logs/${sample_base}_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        2>&1 | tee logs/${sample_base}_picard.log
    
    samtools index $dedup_bam
    
    # Step 5: Generate alignment statistics
    echo "[$(date)] Generating alignment statistics for $sample_base..."
    samtools flagstat $dedup_bam > qc_reports/${sample_base}_flagstat.txt
    samtools stats $dedup_bam > qc_reports/${sample_base}_stats.txt
    
    # Step 6: Create BigWig file using deepTools
    echo "[$(date)] Creating BigWig file for $sample_base..."
    bamCoverage -b $dedup_bam \
        -o $bigwig \
        --binSize 10 \
        --normalizeUsing RPKM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS_PER_JOB \
        2>&1 | tee logs/${sample_base}_bamCoverage.log
    
    # Step 7: Optional - Create normalized BigWig (CPM)
    echo "[$(date)] Creating CPM normalized BigWig for $sample_base..."
    bamCoverage -b $dedup_bam \
        -o ${bigwig/.bw/.CPM.bw} \
        --binSize 10 \
        --normalizeUsing CPM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS_PER_JOB
    
    echo "[$(date)] Completed processing $sample_base"
    echo ""
}

# Export the function and variables so they're available to parallel processes
export -f process_sample
export THREADS_PER_JOB GENOME_INDEX CHROM_SIZES

# Function to download hg38 chromosome sizes if not available
get_chrom_sizes() {
    if [ ! -f "$CHROM_SIZES" ]; then
        echo "Downloading hg38 chromosome sizes..."
        wget -O $CHROM_SIZES https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
    fi
}

# Main processing
main() {
    echo "Cut&Tag to BigWig Processing Pipeline - Parallel Version"
    echo "======================================================"
    echo "Genome index: $GENOME_INDEX"
    echo "Max parallel jobs: $MAX_JOBS"
    echo "Threads per job: $THREADS_PER_JOB"
    echo "Total cores used: $((MAX_JOBS * THREADS_PER_JOB))"
    echo ""
    
    # Check for required tools
    echo "Checking required tools..."
    for tool in bowtie2 samtools picard fastqc bamCoverage; do
        if ! command -v $tool &> /dev/null; then
            echo "ERROR: $tool is not installed or not in PATH"
            exit 1
        fi
    done
    
    # Check for GNU parallel
    if ! command -v parallel &> /dev/null; then
        echo "WARNING: GNU parallel not found. Installing or using alternative method..."
        echo "You can install it with: sudo apt-get install parallel (Ubuntu/Debian)"
        echo "or: brew install parallel (macOS)"
        echo ""
        echo "Falling back to background processes method..."
        USE_PARALLEL=false
    else
        echo "GNU parallel found."
        USE_PARALLEL=true
    fi
    
    echo "All required tools found."
    echo ""
    
    # Get chromosome sizes if needed
    get_chrom_sizes
    
    # Find all R1 files
    r1_files=(raw_data/*_R1.fastq.gz)
    
    if [ ${#r1_files[@]} -eq 0 ]; then
        echo "ERROR: No R1 fastq.gz files found in raw_data/"
        exit 1
    fi
    
    echo "Found ${#r1_files[@]} samples to process"
    echo ""
    
    # Process samples in parallel
    if [ "$USE_PARALLEL" = true ]; then
        echo "Using GNU parallel to process samples..."
        printf '%s\n' "${r1_files[@]}" | parallel -j $MAX_JOBS process_sample
    else
        echo "Using background processes to process samples..."
        
        # Counter for active jobs
        active_jobs=0
        
        for r1 in "${r1_files[@]}"; do
            # Wait if we've reached the maximum number of parallel jobs
            while [ $active_jobs -ge $MAX_JOBS ]; do
                sleep 5
                # Count running background jobs
                active_jobs=$(jobs -r | wc -l)
            done
            
            # Start processing in background
            echo "Starting background process for: $(basename $r1)"
            process_sample "$r1" &
            active_jobs=$((active_jobs + 1))
        done
        
        # Wait for all background jobs to complete
        echo "Waiting for all samples to complete..."
        wait
    fi
    
    echo ""
    echo "All samples processed. Generating reports..."
    
    # Generate MultiQC report
    if command -v multiqc &> /dev/null; then
        echo "Generating MultiQC report..."
        multiqc -f -o qc_reports logs/ qc_reports/
    fi
    
    echo "====================================="
    echo "Pipeline completed!"
    echo "BigWig files are in: bigwig_files/"
    echo "QC reports are in: qc_reports/"
    echo "====================================="
}

# Run the pipeline
main

# Generate a summary table
echo "Generating summary table..."
{
    echo "Sample,Total_Reads,Mapped_Reads,Mapping_Rate,Duplicates,Final_Reads"
    for stats_file in qc_reports/*_flagstat.txt; do
        if [ -f "$stats_file" ]; then
            sample=$(basename $stats_file _flagstat.txt)
            total=$(grep "in total" $stats_file | cut -d' ' -f1)
            mapped=$(grep "mapped (" $stats_file | head -1 | cut -d' ' -f1)
            rate=$(grep "mapped (" $stats_file | head -1 | grep -oP '\(\K[^%]+')
            
            # Get duplicate info from Picard metrics
            dup_file="logs/${sample}_dup_metrics.txt"
            if [ -f "$dup_file" ]; then
                dup_rate=$(grep -A 2 "PERCENT_DUPLICATION" $dup_file | tail -1 | cut -f9)
            else
                dup_rate="N/A"
            fi
            
            echo "$sample,$total,$mapped,$rate%,$dup_rate,$mapped"
        fi
    done
} > qc_reports/processing_summary.csv

echo "Summary saved to: qc_reports/processing_summary.csv"

# Display resource usage summary
echo ""
echo "Resource Usage Summary:"
echo "======================"
echo "Total cores available: 20"
echo "Cores used: $((MAX_JOBS * THREADS_PER_JOB))"
echo "Parallel samples: $MAX_JOBS"
echo "Threads per sample: $THREADS_PER_JOB"
echo "Estimated speedup: ${#r1_files[@]}x faster than sequential processing" 