#!/bin/bash

# Cut&Tag to BigWig Processing Pipeline
# Processes all downloaded Cut&Tag data to create BigWig files

# Configuration
THREADS=8
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
    echo "Processing: $sample_base"
    echo "GSM: $gsm"
    echo "========================================"
    
    # Define file paths
    local r2=${r1/_R1.fastq.gz/_R2.fastq.gz}
    local bam="aligned_bam/${sample_base}.bam"
    local sorted_bam="aligned_bam/${sample_base}.sorted.bam"
    local dedup_bam="aligned_bam/${sample_base}.dedup.bam"
    local bigwig="bigwig_files/${sample_base}.bw"
    
    # Step 1: FastQC on raw reads
    echo "[$(date)] Running FastQC..."
    fastqc -t $THREADS -o qc_reports $r1 $r2 2>&1 | tee logs/${sample_base}_fastqc.log
    
    # Step 2: Align with Bowtie2
    echo "[$(date)] Aligning with Bowtie2..."
    bowtie2 -x $GENOME_INDEX \
        -1 $r1 \
        -2 $r2 \
        -p $THREADS \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        -I 10 -X 1000 \
        2> logs/${sample_base}_bowtie2.log \
        | samtools view -@ $THREADS -bS - > $bam
    
    # Step 3: Sort BAM
    echo "[$(date)] Sorting BAM..."
    samtools sort -@ $THREADS -o $sorted_bam $bam
    samtools index $sorted_bam
    
    # Remove unsorted BAM to save space
    rm $bam
    
    # Step 4: Mark and remove duplicates
    echo "[$(date)] Removing duplicates..."
    picard MarkDuplicates \
        I=$sorted_bam \
        O=$dedup_bam \
        M=logs/${sample_base}_dup_metrics.txt \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        2>&1 | tee logs/${sample_base}_picard.log
    
    samtools index $dedup_bam
    
    # Step 5: Generate alignment statistics
    echo "[$(date)] Generating alignment statistics..."
    samtools flagstat $dedup_bam > qc_reports/${sample_base}_flagstat.txt
    samtools stats $dedup_bam > qc_reports/${sample_base}_stats.txt
    
    # Step 6: Create BigWig file using deepTools
    echo "[$(date)] Creating BigWig file..."
    bamCoverage -b $dedup_bam \
        -o $bigwig \
        --binSize 10 \
        --normalizeUsing RPKM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS \
        2>&1 | tee logs/${sample_base}_bamCoverage.log
    
    # Step 7: Optional - Create normalized BigWig (CPM)
    echo "[$(date)] Creating CPM normalized BigWig..."
    bamCoverage -b $dedup_bam \
        -o ${bigwig/.bw/.CPM.bw} \
        --binSize 10 \
        --normalizeUsing CPM \
        --extendReads \
        --ignoreDuplicates \
        -p $THREADS
    
    echo "[$(date)] Completed processing $sample_base"
    echo ""
}

# Function to download hg38 chromosome sizes if not available
get_chrom_sizes() {
    if [ ! -f "$CHROM_SIZES" ]; then
        echo "Downloading hg38 chromosome sizes..."
        wget -O $CHROM_SIZES https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
    fi
}

# Main processing
main() {
    echo "Cut&Tag to BigWig Processing Pipeline"
    echo "====================================="
    echo "Genome index: $GENOME_INDEX"
    echo "Threads: $THREADS"
    echo ""
    
    # Check for required tools
    echo "Checking required tools..."
    for tool in bowtie2 samtools picard fastqc bamCoverage; do
        if ! command -v $tool &> /dev/null; then
            echo "ERROR: $tool is not installed or not in PATH"
            exit 1
        fi
    done
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
    
    # Process each sample
    for r1 in "${r1_files[@]}"; do
        process_sample "$r1"
    done
    
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

# Optional: Generate a summary table
echo "Generating summary table..."
{
    echo "Sample,Total_Reads,Mapped_Reads,Mapping_Rate,Duplicates,Final_Reads"
    for stats_file in qc_reports/*_flagstat.txt; do
        sample=$(basename $stats_file _flagstat.txt)
        total=$(grep "in total" $stats_file | cut -d' ' -f1)
        mapped=$(grep "mapped (" $stats_file | cut -d' ' -f1)
        rate=$(grep "mapped (" $stats_file | grep -oP '\(\K[^%]+')
        
        # Get duplicate info from Picard metrics
        dup_file="logs/${sample}_dup_metrics.txt"
        if [ -f "$dup_file" ]; then
            dup_rate=$(grep -A 2 "PERCENT_DUPLICATION" $dup_file | tail -1 | cut -f9)
        else
            dup_rate="N/A"
        fi
        
        echo "$sample,$total,$mapped,$rate%,$dup_rate,$mapped"
    done
} > qc_reports/processing_summary.csv

echo "Summary saved to: qc_reports/processing_summary.csv"