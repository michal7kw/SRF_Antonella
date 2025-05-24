#!/bin/bash

# Cut&Tag to BigWig Processing Pipeline - Fixed Parallel Version
# Processes all downloaded Cut&Tag data to create BigWig files using parallel processing

# Configuration
MAX_JOBS=5          # Number of samples to process in parallel (leave some cores for system)
THREADS_PER_JOB=4   # Threads per sample (5 jobs * 4 threads = 20 cores)
GENOME_INDEX="/home/michal/COMMON_FILES/bowtie2_indexes/hg38/GRCh38_noalt_as/GRCh38_noalt_as" 
CHROM_SIZES="/home/michal/COMMON_FILES/hg38.chrom.sizes"

# Create output directories
mkdir -p aligned_bam
mkdir -p bigwig_files
mkdir -p qc_reports
mkdir -p logs

# Function to check Java version and suggest alternatives
check_java_version() {
    echo "Checking Java version for Picard compatibility..."
    
    # If variables are already set (e.g., by wrapper script), use them
    if [[ -n "${USE_SAMTOOLS_MARKDUP:-}" ]]; then
        echo "Using pre-configured duplicate removal method:"
        if [[ "$USE_SAMTOOLS_MARKDUP" == "true" ]]; then
            echo "✅ samtools markdup (Java version compatible)"
        else
            echo "✅ Picard MarkDuplicates (Java 17+ required)"
        fi
        return 0
    fi
    
    if command -v java &> /dev/null; then
        JAVA_VERSION=$(java -version 2>&1 | grep -oP 'version "([0-9]+)' | grep -oP '[0-9]+' | head -1)
        echo "Java version detected: $JAVA_VERSION"
        
        if [ "$JAVA_VERSION" -lt 17 ]; then
            echo "WARNING: Picard requires Java 17+, but Java $JAVA_VERSION detected."
            echo "Options:"
            echo "1. Use samtools markdup instead of Picard (recommended)"
            echo "2. Install Java 17+ (conda install openjdk=17)"
            echo "3. Skip duplicate removal (not recommended)"
            echo ""
            read -p "Choose option (1/2/3): " CHOICE
            case $CHOICE in
                1) USE_SAMTOOLS_MARKDUP=true ;;
                2) echo "Please install Java 17+ and re-run the script"; exit 1 ;;
                3) SKIP_DEDUP=true ;;
                *) echo "Invalid choice, using samtools markdup"; USE_SAMTOOLS_MARKDUP=true ;;
            esac
        else
            USE_SAMTOOLS_MARKDUP=false
        fi
    else
        echo "Java not found! Installing via conda..."
        conda install -y openjdk=17
        USE_SAMTOOLS_MARKDUP=false
    fi
}

# Function to process a single sample with error handling
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
    
    # Check if input files exist
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "ERROR: Input files not found for $sample_base"
        echo "R1: $r1"
        echo "R2: $r2"
        return 1
    fi
    
    # Step 1: FastQC on raw reads
    echo "[$(date)] Running FastQC for $sample_base..."
    if ! fastqc -t $THREADS_PER_JOB -o qc_reports $r1 $r2 2>&1 | tee logs/${sample_base}_fastqc.log; then
        echo "ERROR: FastQC failed for $sample_base"
        return 1
    fi
    
    # Step 2: Align with Bowtie2
    echo "[$(date)] Aligning $sample_base with Bowtie2..."
    
    # Check if genome index exists
    if [[ ! -f "${GENOME_INDEX}.1.bt2" ]]; then
        echo "ERROR: Bowtie2 index not found at $GENOME_INDEX"
        echo "Please check the path or create the index"
        return 1
    fi
    
    if ! bowtie2 -x $GENOME_INDEX \
        -1 $r1 \
        -2 $r2 \
        -p $THREADS_PER_JOB \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -I 10 -X 1000 \
        2> logs/${sample_base}_bowtie2.log \
        | samtools view -@ $THREADS_PER_JOB -bS - > $bam; then
        echo "ERROR: Bowtie2 alignment failed for $sample_base"
        cat logs/${sample_base}_bowtie2.log
        return 1
    fi
    
    # Check if BAM file was created
    if [[ ! -f "$bam" || ! -s "$bam" ]]; then
        echo "ERROR: BAM file not created or empty for $sample_base"
        return 1
    fi
    
    # Step 3: Sort BAM
    echo "[$(date)] Sorting BAM for $sample_base..."
    if ! samtools sort -@ $THREADS_PER_JOB -o $sorted_bam $bam; then
        echo "ERROR: BAM sorting failed for $sample_base"
        return 1
    fi
    
    if ! samtools index $sorted_bam; then
        echo "ERROR: BAM indexing failed for $sample_base"
        return 1
    fi
    
    # Remove unsorted BAM to save space
    rm $bam
    
    # Step 4: Mark and remove duplicates
    echo "[$(date)] Removing duplicates for $sample_base..."
    
    if [[ "$USE_SAMTOOLS_MARKDUP" == "true" ]]; then
        echo "Using samtools markdup (Java version compatible)..."
        
        # samtools markdup workflow: sort by name → fixmate → sort by coordinate → markdup
        local name_sorted="aligned_bam/${sample_base}.name_sorted.bam"
        local fixmate_bam="aligned_bam/${sample_base}.fixmate.bam"
        local fixmate_sorted="aligned_bam/${sample_base}.fixmate.sorted.bam"
        
        echo "Sorting BAM by query name for fixmate..."
        if ! samtools sort -n -@ $THREADS_PER_JOB -o $name_sorted $sorted_bam; then
            echo "ERROR: name sorting failed for $sample_base"
            echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
            cp $sorted_bam $dedup_bam
        else
            echo "Running samtools fixmate..."
            if ! samtools fixmate -m $name_sorted $fixmate_bam; then
                echo "ERROR: samtools fixmate failed for $sample_base"
                echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                cp $sorted_bam $dedup_bam
                rm -f $name_sorted
            else
                # Sort the fixmate BAM by coordinate (required for markdup)
                echo "Sorting fixmate BAM by coordinate..."
                if ! samtools sort -@ $THREADS_PER_JOB -o $fixmate_sorted $fixmate_bam; then
                    echo "ERROR: fixmate BAM sorting failed for $sample_base"
                    echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                    cp $sorted_bam $dedup_bam
                    rm -f $name_sorted $fixmate_bam
                else
                    # Now run markdup
                    echo "Running samtools markdup..."
                    if ! samtools markdup -r -@ $THREADS_PER_JOB -s $fixmate_sorted $dedup_bam 2>&1 | tee logs/${sample_base}_markdup.log; then
                        echo "ERROR: samtools markdup failed for $sample_base"
                        echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                        cp $sorted_bam $dedup_bam
                    fi
                    
                    # Clean up intermediate files
                    rm -f $name_sorted $fixmate_bam $fixmate_sorted
                fi
            fi
        fi
    elif [[ "$SKIP_DEDUP" == "true" ]]; then
        echo "Skipping duplicate removal..."
        cp $sorted_bam $dedup_bam
    else
        echo "Using Picard MarkDuplicates..."
        if ! picard MarkDuplicates \
            I=$sorted_bam \
            O=$dedup_bam \
            M=logs/${sample_base}_dup_metrics.txt \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tee logs/${sample_base}_picard.log; then
            echo "ERROR: Picard MarkDuplicates failed for $sample_base"
            echo "Falling back to samtools markdup..."
            
            # Fallback with proper fixmate workflow
            local name_sorted="aligned_bam/${sample_base}.name_sorted.bam"
            local fixmate_bam="aligned_bam/${sample_base}.fixmate.bam"
            local fixmate_sorted="aligned_bam/${sample_base}.fixmate.sorted.bam"
            
            echo "Sorting BAM by query name for fixmate (fallback)..."
            if ! samtools sort -n -@ $THREADS_PER_JOB -o $name_sorted $sorted_bam; then
                echo "ERROR: name sorting fallback failed for $sample_base"
                echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                cp $sorted_bam $dedup_bam
            else
                echo "Running samtools fixmate (fallback)..."
                if ! samtools fixmate -m $name_sorted $fixmate_bam; then
                    echo "ERROR: samtools fixmate fallback failed for $sample_base"
                    echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                    cp $sorted_bam $dedup_bam
                    rm -f $name_sorted
                else
                    echo "Sorting fixmate BAM by coordinate (fallback)..."
                    if ! samtools sort -@ $THREADS_PER_JOB -o $fixmate_sorted $fixmate_bam; then
                        echo "ERROR: fixmate BAM sorting fallback failed for $sample_base"
                        echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                        cp $sorted_bam $dedup_bam
                        rm -f $name_sorted $fixmate_bam
                    else
                        if ! samtools markdup -r -@ $THREADS_PER_JOB -s $fixmate_sorted $dedup_bam 2>&1 | tee logs/${sample_base}_markdup_fallback.log; then
                            echo "ERROR: Fallback markdup also failed for $sample_base"
                            echo "WARNING: Skipping duplicate removal and continuing with sorted BAM"
                            cp $sorted_bam $dedup_bam
                        fi
                        
                        # Clean up intermediate files
                        rm -f $name_sorted $fixmate_bam $fixmate_sorted
                    fi
                fi
            fi
        fi
    fi
    
    if ! samtools index $dedup_bam; then
        echo "ERROR: Dedup BAM indexing failed for $sample_base"
        return 1
    fi
    
    # Check if dedup BAM has any reads
    read_count=$(samtools view -c $dedup_bam)
    if [[ $read_count -eq 0 ]]; then
        echo "ERROR: No reads in deduplicated BAM file for $sample_base"
        echo "This suggests alignment failed or all reads were duplicates"
        return 1
    fi
    
    echo "Deduplicated BAM contains $read_count reads"
    
    # Step 5: Generate alignment statistics
    echo "[$(date)] Generating alignment statistics for $sample_base..."
    samtools flagstat $dedup_bam > qc_reports/${sample_base}_flagstat.txt
    samtools stats $dedup_bam > qc_reports/${sample_base}_stats.txt
    
    # Step 6: Create BigWig file using deepTools
    echo "[$(date)] Creating BigWig file for $sample_base..."
    if ! bamCoverage -b $dedup_bam \
        -o $bigwig \
        --binSize 10 \
        --normalizeUsing RPKM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS_PER_JOB \
        2>&1 | tee logs/${sample_base}_bamCoverage.log; then
        echo "ERROR: BigWig creation failed for $sample_base"
        return 1
    fi
    
    # Step 7: Create normalized BigWig (CPM)
    echo "[$(date)] Creating CPM normalized BigWig for $sample_base..."
    if ! bamCoverage -b $dedup_bam \
        -o ${bigwig/.bw/.CPM.bw} \
        --binSize 10 \
        --normalizeUsing CPM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS_PER_JOB; then
        echo "WARNING: CPM BigWig creation failed for $sample_base (continuing...)"
    fi
    
    echo "[$(date)] Successfully completed processing $sample_base"
    echo ""
    return 0
}

# Export the function and variables so they're available to parallel processes
export -f process_sample
export THREADS_PER_JOB GENOME_INDEX CHROM_SIZES USE_SAMTOOLS_MARKDUP SKIP_DEDUP

# Function to download hg38 chromosome sizes if not available
get_chrom_sizes() {
    if [ ! -f "$CHROM_SIZES" ]; then
        echo "Downloading hg38 chromosome sizes..."
        wget -O $CHROM_SIZES https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
    fi
}

# Main processing
main() {
    echo "Cut&Tag to BigWig Processing Pipeline - Fixed Parallel Version"
    echo "============================================================="
    echo "Genome index: $GENOME_INDEX"
    echo "Max parallel jobs: $MAX_JOBS"
    echo "Threads per job: $THREADS_PER_JOB"
    echo "Total cores used: $((MAX_JOBS * THREADS_PER_JOB))"
    echo ""
    
    # Check Java version and set duplicate removal method
    check_java_version
    
    # Enable strict error handling after Java version check
    set -euo pipefail
    
    # Check for required tools
    echo "Checking required tools..."
    MISSING_TOOLS=()
    
    for tool in bowtie2 samtools fastqc bamCoverage; do
        if ! command -v $tool &> /dev/null; then
            MISSING_TOOLS+=($tool)
        fi
    done
    
    if [ ${#MISSING_TOOLS[@]} -ne 0 ]; then
        echo "ERROR: Missing required tools: ${MISSING_TOOLS[*]}"
        echo "Please install them using conda:"
        echo "conda install -c bioconda ${MISSING_TOOLS[*]}"
        exit 1
    fi
    
    echo "All required tools found."
    echo ""
    
    # Check if genome index exists
    if [[ ! -f "${GENOME_INDEX}.1.bt2" ]]; then
        echo "ERROR: Bowtie2 index not found at $GENOME_INDEX"
        echo "Available options:"
        echo "1. Create index: bowtie2-build genome.fa $GENOME_INDEX"
        echo "2. Download pre-built index"
        echo "3. Update GENOME_INDEX path in script"
        exit 1
    fi
    
    # Get chromosome sizes if needed
    get_chrom_sizes
    
    # Find all R1 files
    r1_files=(raw_data/*_R1.fastq.gz)
    
    if [ ${#r1_files[@]} -eq 0 ]; then
        echo "ERROR: No R1 fastq.gz files found in raw_data/"
        exit 1
    fi
    
    echo "Found ${#r1_files[@]} samples to process:"
    for r1 in "${r1_files[@]}"; do
        echo "  - $(basename $r1)"
    done
    echo ""
    
    # Process samples in parallel with error tracking
    echo "Starting parallel processing..."
    FAILED_SAMPLES=()
    
    # Check for GNU parallel
    if command -v parallel &> /dev/null; then
        echo "Using GNU parallel to process samples..."
        
        # Create temporary file to track failures
        FAILED_FILE=$(mktemp)
        
        printf '%s\n' "${r1_files[@]}" | parallel -j $MAX_JOBS --halt now,fail=1 "process_sample {} || echo {} >> $FAILED_FILE"
        
        # Check for failures
        if [[ -f "$FAILED_FILE" && -s "$FAILED_FILE" ]]; then
            echo "ERROR: Some samples failed processing:"
            cat "$FAILED_FILE"
            rm "$FAILED_FILE"
            exit 1
        fi
        rm -f "$FAILED_FILE"
        
    else
        echo "Using background processes to process samples..."
        
        # Process samples with error tracking
        for r1 in "${r1_files[@]}"; do
            # Wait if we've reached the maximum number of parallel jobs
            while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
                sleep 5
            done
            
            echo "Starting background process for: $(basename $r1)"
            {
                if ! process_sample "$r1"; then
                    echo "FAILED: $r1" >> processing_failures.txt
                fi
            } &
        done
        
        # Wait for all background jobs to complete
        echo "Waiting for all samples to complete..."
        wait
        
        # Check for failures
        if [[ -f "processing_failures.txt" ]]; then
            echo "ERROR: Some samples failed processing:"
            cat processing_failures.txt
            exit 1
        fi
    fi
    
    echo ""
    echo "All samples processed successfully!"
    
    # Generate MultiQC report
    if command -v multiqc &> /dev/null; then
        echo "Generating MultiQC report..."
        multiqc -f -o qc_reports logs/ qc_reports/
    fi
    
    echo "====================================="
    echo "Pipeline completed successfully!"
    echo "BigWig files are in: bigwig_files/"
    echo "QC reports are in: qc_reports/"
    echo "====================================="
}

# Run the pipeline
main "$@"

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
            
            # Get duplicate info
            if [[ "$USE_SAMTOOLS_MARKDUP" == "true" ]]; then
                dup_rate="N/A (samtools)"
            else
                dup_file="logs/${sample}_dup_metrics.txt"
                if [ -f "$dup_file" ]; then
                    dup_rate=$(grep -A 2 "PERCENT_DUPLICATION" $dup_file | tail -1 | cut -f9)
                else
                    dup_rate="N/A"
                fi
            fi
            
            echo "$sample,$total,$mapped,$rate%,$dup_rate,$mapped"
        fi
    done
} > qc_reports/processing_summary.csv

echo "Summary saved to: qc_reports/processing_summary.csv"

# Display final summary
echo ""
echo "Final Summary:"
echo "============="
echo "✓ Samples processed: $(ls bigwig_files/*.bw 2>/dev/null | wc -l)"
echo "✓ BigWig files created: $(ls bigwig_files/ 2>/dev/null | wc -l)"
echo "✓ QC reports: $(ls qc_reports/ 2>/dev/null | wc -l)"
echo "✓ Log files: $(ls logs/ 2>/dev/null | wc -l)" 