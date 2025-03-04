#!/bin/bash
#SBATCH --job-name=create_bigwig_cutntag
#SBATCH --output=logs/create_bigwig_cutntag_%a.out
#SBATCH --error=logs/create_bigwig_cutntag_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --array=0-5

# Script to create BigWig files from FASTQ files for SRF_SES_V5 project
# These files are from a study on glioblastoma cells (SNB19) with SES and SOX2 transduction
# Optimized for Cut&Tag data analysis
# Modified to reuse already correctly computed results

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/data_from_ncbi_corrected

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create necessary directories
mkdir -p logs
mkdir -p bam
mkdir -p bigwig
mkdir -p peaks

# Reference genome path - updated with correct path
REFERENCE_GENOME="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# Bowtie2 index path - updated with correct path (without file extension)
BOWTIE2_INDEX="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/GRCh38"
# Chromosome sizes file for bigWig generation
CHROM_SIZES="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes"

# Array of SRR accessions to process
declare -a SRR_ACCESSIONS=(
    "SRR18590288"  # SES-transduced SNB19
    "SRR18590287"  # SES-transduced SNB19
    "SRR18590286"  # SES-transduced SNB19
    "SRR18590285"  # SOX2-transduced SNB19
    "SRR18590284"  # SOX2-transduced SNB19
    "SRR18590283"  # SOX2-transduced SNB19
)

# Get the current SRR accession based on the array task ID
current_accession=${SRR_ACCESSIONS[$SLURM_ARRAY_TASK_ID]}

echo "Starting Cut&Tag analysis for: $current_accession"

# Check if the FASTQ files exist
if [ ! -f "./fastq/${current_accession}_1.fastq.gz" ] || [ ! -f "./fastq/${current_accession}_2.fastq.gz" ]; then
    echo "FASTQ files for $current_accession not found. Please run sra_to_fastq.sh first."
    exit 1
fi

# Get sample information from metadata
treatment=$(grep "$current_accession" sample_metadata.tsv | cut -f5)
if [[ "$treatment" == "SES-transduced" ]]; then
    # Find which replicate this is among SES samples
    ses_samples=($(grep "SES-transduced" sample_metadata.tsv | cut -f1))
    for i in "${!ses_samples[@]}"; do
        if [[ "${ses_samples[$i]}" == "$current_accession" ]]; then
            replicate=$((i+1))
            break
        fi
    done
    sample_name="SES_rep${replicate}"
elif [[ "$treatment" == "SOX2-transduced" ]]; then
    # Find which replicate this is among SOX2 samples
    sox2_samples=($(grep "SOX2-transduced" sample_metadata.tsv | cut -f1))
    for i in "${!sox2_samples[@]}"; do
        if [[ "${sox2_samples[$i]}" == "$current_accession" ]]; then
            replicate=$((i+1))
            break
        fi
    done
    sample_name="SOX2_rep${replicate}"
else
    sample_name="${current_accession}"
fi

echo "Processing sample: $sample_name (SRR: $current_accession)"

# Function to check file age in minutes
file_age_minutes() {
    local file=$1
    local current_time=$(date +%s)
    local file_time=$(stat -c %Y "$file")
    local age_seconds=$((current_time - file_time))
    echo $((age_seconds / 60))
}

# Function to check if a file exists and is not empty
file_exists_and_not_empty() {
    local file=$1
    if [ -f "$file" ] && [ -s "$file" ]; then
        return 0  # True
    else
        return 1  # False
    fi
}

# Step 1: Trim adapters with Trimmomatic (important for Cut&Tag data)
echo "Checking for existing trimmed FASTQ files..."

# Check if trimmed files already exist and are not empty
if file_exists_and_not_empty "./fastq/${sample_name}_R1_trimmed_paired.fastq.gz" && \
   file_exists_and_not_empty "./fastq/${sample_name}_R2_trimmed_paired.fastq.gz"; then
    echo "Trimmed FASTQ files already exist for $sample_name. Skipping trimming step."
else
    echo "Trimming adapters with Trimmomatic..."

    # Check if Trimmomatic is installed
    if ! command -v trimmomatic &> /dev/null; then
        echo "Trimmomatic is not installed or not in PATH. Please install it first."
        exit 1
    fi

    # Path to Trimmomatic adapter file
    ADAPTER_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/TruSeq3-PE.fa"

    trimmomatic PE -threads 8 \
        ./fastq/${current_accession}_1.fastq.gz \
        ./fastq/${current_accession}_2.fastq.gz \
        ./fastq/${sample_name}_R1_trimmed_paired.fastq.gz \
        ./fastq/${sample_name}_R1_trimmed_unpaired.fastq.gz \
        ./fastq/${sample_name}_R2_trimmed_paired.fastq.gz \
        ./fastq/${sample_name}_R2_trimmed_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    if [ $? -ne 0 ]; then
        echo "Error during adapter trimming with Trimmomatic."
        exit 1
    fi
fi

# Step 2: Align reads to reference genome using Bowtie2
echo "Checking for existing alignment files..."

# Check if sorted BAM file already exists and is not empty
if file_exists_and_not_empty "./bam/${sample_name}.sorted.bam" && \
   file_exists_and_not_empty "./bam/${sample_name}.sorted.bam.bai"; then
    echo "Sorted BAM file already exists for $sample_name. Skipping alignment step."
else
    echo "Aligning reads to reference genome using Bowtie2..."

    # Check if Bowtie2 is installed
    if ! command -v bowtie2 &> /dev/null; then
        echo "Bowtie2 is not installed or not in PATH. Please install it first."
        exit 1
    fi

    # Check if samtools is installed
    if ! command -v samtools &> /dev/null; then
        echo "samtools is not installed or not in PATH. Please install it first."
        exit 1
    fi

    # Align with parameters optimized for Cut&Tag data
    bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        --phred33 -I 10 -X 700 \
        -p 8 \
        -x ${BOWTIE2_INDEX} \
        -1 ./fastq/${sample_name}_R1_trimmed_paired.fastq.gz \
        -2 ./fastq/${sample_name}_R2_trimmed_paired.fastq.gz \
        -S ./bam/${sample_name}.sam \
        2> ./bam/${sample_name}_alignment_summary.txt

    if [ $? -ne 0 ]; then
        echo "Error during alignment with Bowtie2."
        exit 1
    fi

    # Convert to BAM
    echo "Converting SAM to BAM, sorting, and filtering..."
    samtools view -@ 8 -bS ./bam/${sample_name}.sam > ./bam/${sample_name}.bam

    # Filter for high-quality alignments and proper pairs (important for Cut&Tag)
    samtools view -@ 8 -b -f 2 -q 10 ./bam/${sample_name}.bam > ./bam/${sample_name}.filtered.bam

    # Sort and index
    samtools sort -@ 8 ./bam/${sample_name}.filtered.bam -o ./bam/${sample_name}.sorted.bam
    samtools index ./bam/${sample_name}.sorted.bam

    # Remove intermediate files to save space
    rm ./bam/${sample_name}.sam
    rm ./bam/${sample_name}.bam
    rm ./bam/${sample_name}.filtered.bam
fi

# Step 3: Mark duplicates
echo "Checking for existing deduplicated BAM files..."

# Check if deduplicated BAM file already exists and is not empty
if file_exists_and_not_empty "./bam/${sample_name}.dedup.bam" && \
   file_exists_and_not_empty "./bam/${sample_name}.dedup.bam.bai"; then
    echo "Deduplicated BAM file already exists for $sample_name. Skipping deduplication step."
else
    echo "Marking duplicates with samtools..."

    # Create a name-sorted BAM file (required for markdup)
    samtools sort -@ 8 -n ./bam/${sample_name}.sorted.bam -o ./bam/${sample_name}.namesorted.bam

    # Mark duplicates with samtools
    samtools fixmate -@ 8 -m ./bam/${sample_name}.namesorted.bam ./bam/${sample_name}.fixmate.bam
    samtools sort -@ 8 ./bam/${sample_name}.fixmate.bam -o ./bam/${sample_name}.positionsorted.bam
    samtools markdup -@ 8 -r ./bam/${sample_name}.positionsorted.bam ./bam/${sample_name}.dedup.bam

    # Index the final BAM file
    samtools index ./bam/${sample_name}.dedup.bam

    # Clean up intermediate files
    rm ./bam/${sample_name}.namesorted.bam
    rm ./bam/${sample_name}.fixmate.bam
    rm ./bam/${sample_name}.positionsorted.bam

    # Verify the deduped BAM file exists before proceeding
    if [ ! -f "./bam/${sample_name}.dedup.bam" ]; then
        echo "Error: Deduplication failed. The deduped BAM file does not exist."
        exit 1
    fi
fi

# Step 4: Generate BigWig files
echo "Checking for existing BigWig files..."

# Check if BigWig file already exists and is not empty
if file_exists_and_not_empty "./bigwig/${sample_name}.bw"; then
    echo "BigWig file already exists for $sample_name. Skipping BigWig generation step."
else
    echo "Generating BigWig files..."

    # Check if deepTools is installed
    if ! command -v bamCoverage &> /dev/null; then
        echo "deepTools is not installed or not in PATH. Please install it first."
        exit 1
    fi

    # Generate normalized BigWig file (RPGC normalization is better for Cut&Tag)
    # Also use a smaller bin size for higher resolution
    bamCoverage --bam ./bam/${sample_name}.dedup.bam \
        --outFileName ./bigwig/${sample_name}.bw \
        --outFileFormat bigwig \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --binSize 5 \
        --smoothLength 15 \
        --numberOfProcessors 8 \
        --ignoreForNormalization chrX chrY chrM \
        --extendReads

    # If bamCoverage fails, try an alternative approach using bedtools and bedGraphToBigWig
    if [ $? -ne 0 ]; then
        echo "bamCoverage failed. Trying alternative approach with bedtools and bedGraphToBigWig..."
        
        # Check if bedtools is installed
        if ! command -v bedtools &> /dev/null; then
            echo "bedtools is not installed or not in PATH. Cannot proceed with alternative approach."
            exit 1
        fi
        
        # Check if bedGraphToBigWig is installed
        if ! command -v bedGraphToBigWig &> /dev/null; then
            echo "bedGraphToBigWig is not installed or not in PATH. Cannot proceed with alternative approach."
            exit 1
        fi
        
        # Generate bedGraph
        bedtools genomecov -ibam ./bam/${sample_name}.dedup.bam -bg -scale 1.0 > ./bigwig/${sample_name}.bedgraph
        
        # Sort bedGraph (required for bedGraphToBigWig)
        sort -k1,1 -k2,2n ./bigwig/${sample_name}.bedgraph > ./bigwig/${sample_name}.sorted.bedgraph
        
        # Convert to bigWig
        bedGraphToBigWig ./bigwig/${sample_name}.sorted.bedgraph ${CHROM_SIZES} ./bigwig/${sample_name}.bw
        
        # Clean up
        rm ./bigwig/${sample_name}.bedgraph
        rm ./bigwig/${sample_name}.sorted.bedgraph
        
        if [ -f "./bigwig/${sample_name}.bw" ]; then
            echo "Successfully created BigWig file for $sample_name using alternative approach"
        else
            echo "Error creating BigWig file for $sample_name using alternative approach"
            exit 1
        fi
    else
        echo "Successfully created BigWig file for $sample_name"
    fi
fi

# Step 5: Call peaks with MACS2
echo "Checking for existing peak files..."

# Check if peak file already exists and is not empty
if file_exists_and_not_empty "./peaks/${sample_name}_peaks.narrowPeak"; then
    echo "Peak file already exists for $sample_name. Skipping peak calling step."
else
    echo "Calling peaks with MACS2..."

    # Check if MACS2 is installed
    if ! command -v macs2 &> /dev/null; then
        echo "MACS2 is not installed or not in PATH. Please install it first."
        exit 1
    fi

    # Verify the deduped BAM file exists before calling peaks
    if [ ! -f "./bam/${sample_name}.dedup.bam" ]; then
        echo "Error: The deduped BAM file does not exist. Cannot call peaks."
        exit 1
    fi

    # Call peaks with parameters optimized for Cut&Tag data
    # Use --format BAMPE for paired-end data
    # Use a more stringent q-value for Cut&Tag
    macs2 callpeak \
        -t ./bam/${sample_name}.dedup.bam \
        -f BAMPE \
        -g hs \
        -n ${sample_name} \
        --outdir ./peaks \
        -q 0.01 \
        --keep-dup all \
        --call-summits

    if [ $? -eq 0 ]; then
        echo "Successfully called peaks for $sample_name"
    else
        echo "Error calling peaks for $sample_name"
        # Continue execution even if peak calling fails
    fi
fi

# Create or update track line files for IGV
# Check if this sample is already in the track lines file
if ! grep -q "${sample_name}" ./bigwig/track_lines.txt; then
    echo "Adding BigWig track line for $sample_name to track_lines.txt"
    echo "track type=bigWig name=\"$sample_name\" description=\"$treatment replicate $replicate Cut&Tag\" bigDataUrl=./bigwig/${sample_name}.bw" >> ./bigwig/track_lines.txt
fi

# Also create a track for the peaks if they exist and not already in the track lines file
if file_exists_and_not_empty "./peaks/${sample_name}_peaks.narrowPeak"; then
    if ! grep -q "${sample_name}_peaks" ./peaks/track_lines.txt 2>/dev/null; then
        echo "Adding peak track line for $sample_name to peaks/track_lines.txt"
        echo "track type=bed name=\"${sample_name}_peaks\" description=\"$treatment replicate $replicate Cut&Tag peaks\" itemRgb=On" >> ./peaks/track_lines.txt
        
        # Add browser position line only once at the beginning of the file
        if [ ! -s ./peaks/track_lines.txt ] || ! grep -q "browser position" ./peaks/track_lines.txt; then
            sed -i '1i browser position chr1:1-10000' ./peaks/track_lines.txt
        fi
    fi
fi

echo "Cut&Tag analysis completed for $sample_name." 