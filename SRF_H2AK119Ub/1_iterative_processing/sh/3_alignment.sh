#!/bin/bash
#SBATCH --job-name=3_alignment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5  # Process 6 samples in parallel (GFP_1-3 and YAF_1-3)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3_alignment_%a.err"
#SBATCH --output="logs/3_alignment_%a.out"

# This script performs alignment of paired-end CUT&Tag sequencing data to a reference genome
# It processes multiple samples in parallel using a SLURM array job
# For each sample, it:
# 1. Aligns trimmed reads to the genome using Bowtie2 with CUT&Tag-specific parameters
# 2. Filters aligned reads based on mapping quality
# 3. Marks and optionally removes PCR duplicates
# 4. Generates comprehensive QC metrics and visualizations

# Input files (per sample):
# - analysis/trimmed/{sample}_R1_paired.fastq.gz - Trimmed forward reads
# - analysis/trimmed/{sample}_R2_paired.fastq.gz - Trimmed reverse reads
# - analysis/qc/{sample}_trimming_metrics.json - QC metrics from trimming step
# - genome/GRCh38.* - Bowtie2 genome index files

# Output files (per sample):
# - analysis/aligned/{sample}.dedup.bam - Final filtered and deduplicated BAM file
# - analysis/aligned/{sample}.dedup.bam.bai - Index for final BAM
# - analysis/qc/{sample}_mapping_stats.csv - Basic mapping statistics
# - analysis/qc/{sample}_prefilter_flagstat.txt - Pre-filtering alignment stats
# - analysis/qc/{sample}_mapq_distribution.txt - Distribution of mapping qualities
# - analysis/qc/{sample}_prededup_flagstat.txt - Pre-deduplication stats
# - analysis/qc/{sample}_dup_metrics.txt - PCR duplicate metrics
# - analysis/qc/{sample}_alignment_summary.txt - Comprehensive alignment report
# - analysis/qc/{sample}_flagstat.txt - Final alignment statistics
# - analysis/qc/{sample}_insert_size_metrics.txt - Fragment size metrics
# - analysis/qc/{sample}_insert_size_histogram.pdf - Fragment size distribution plot
# - logs/aligned/{sample}.align.log - Bowtie2 alignment log

# Exit script on error, undefined variables, and pipe failures
set -e
set -u
set -o pipefail

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check if a command exists
# Used to verify required tools are available
check_command() {
    if ! command -v $1 &> /dev/null; then
        log_message "ERROR: $1 command not found. Please ensure it is installed and in PATH"
        exit 1
    fi
}

# Check required commands
log_message "Checking required commands..."
check_command bowtie2
check_command samtools
check_command picard
check_command java

# Print conda environment info for reproducibility
log_message "Conda environment info:"
conda info
log_message "Checking required packages:"
for pkg in bowtie2 samtools picard; do
    conda list $pkg
done

# Define working directory and create necessary subdirectories
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
log_message "Changing to working directory: $WORKDIR"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories for output files
log_message "Creating output directories..."
mkdir -p analysis/{aligned,qc} logs/aligned || { log_message "ERROR: Failed to create directories"; exit 1; }

# Function to check if output exists and is not empty
# Used to verify successful generation of output files
check_output() {
    local file=$1
    if [[ ! -s $file ]]; then
        log_message "ERROR: Output file $file is empty or does not exist"
        exit 1
    else
        log_message "Successfully created: $file ($(du -h $file | cut -f1))"
    fi
}

# Function to get mapping statistics
# Calculates and saves basic mapping metrics
get_mapping_stats() {
    local bam=$1
    local prefix=$2
    
    log_message "Calculating mapping statistics for $bam..."
    # Get total reads
    total_reads=$(samtools view -c $bam)
    # Get mapped reads
    mapped_reads=$(samtools view -c -F 4 $bam)
    # Calculate mapping rate
    mapping_rate=$(echo "scale=4; $mapped_reads / $total_reads * 100" | bc)
    
    log_message "Total reads: $total_reads"
    log_message "Mapped reads: $mapped_reads"
    log_message "Mapping rate: ${mapping_rate}%"
    
    echo "${prefix},${total_reads},${mapped_reads},${mapping_rate}" > analysis/qc/${sample}_mapping_stats.csv
}

# Define samples array and get current sample based on array task ID
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

log_message "Processing sample: ${sample}"

# Check if input files exist
if [[ ! -f analysis/trimmed/${sample}_R1_paired.fastq.gz ]] || [[ ! -f analysis/trimmed/${sample}_R2_paired.fastq.gz ]]; then
    log_message "ERROR: Input fastq files not found for ${sample}"
    log_message "Expected files:"
    log_message "  - analysis/trimmed/${sample}_R1_paired.fastq.gz"
    log_message "  - analysis/trimmed/${sample}_R2_paired.fastq.gz"
    exit 1
fi

# Check if trimming metrics exist
if [[ ! -f analysis/qc/${sample}_trimming_metrics.json ]]; then
    log_message "ERROR: Trimming metrics not found. Please run 2a_quality_control.sh first"
    exit 1
fi

# Define genome directory and check index
# GENOME_DIR="${WORKDIR}/genome"
# GENOME_INDEX="${GENOME_DIR}/GRCh38"

GENOME_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/genome"
GENOME_INDEX="${GENOME_DIR}/GRCh38"

if [[ ! -d "${GENOME_DIR}" ]]; then
    log_message "ERROR: Genome directory not found at ${GENOME_DIR}"
    log_message "Please run the genome preparation script first:"
    log_message "sbatch 1_genome_preparation.sh"
    exit 1
fi

if [[ ! -f "${GENOME_INDEX}.1.bt2" ]]; then
    log_message "ERROR: Bowtie2 index not found at ${GENOME_INDEX}"
    log_message "Expected files:"
    log_message "  - ${GENOME_INDEX}.1.bt2"
    log_message "  - ${GENOME_INDEX}.2.bt2"
    log_message "  - ${GENOME_INDEX}.3.bt2"
    log_message "  - ${GENOME_INDEX}.4.bt2"
    log_message "  - ${GENOME_INDEX}.rev.1.bt2"
    log_message "  - ${GENOME_INDEX}.rev.2.bt2"
    log_message "Please run the genome preparation script first:"
    log_message "sbatch 1_genome_preparation.sh"
    exit 1
fi

# Create temporary directory for intermediate files
tmp_dir="analysis/aligned/tmp_${sample}"
log_message "Creating temporary directory: ${tmp_dir}"
mkdir -p ${tmp_dir} || { log_message "ERROR: Failed to create temporary directory"; exit 1; }

# Read trimming metrics for QC reporting
if [[ -f analysis/qc/${sample}_trimming_metrics.json ]]; then
    trim_metrics=$(cat analysis/qc/${sample}_trimming_metrics.json)
    initial_reads=$(echo $trim_metrics | jq -r '.clean_r1_reads')
    log_message "Found trimming metrics: ${initial_reads} input reads"
else
    log_message "WARNING: No trimming metrics found, will calculate statistics from scratch"
fi

# Run Bowtie2 alignment with CUT&Tag-specific parameters:
# - local alignment for better peak detection
# - very-sensitive-local for higher accuracy
# - no mixed/discordant alignments to ensure proper pairs
log_message "Starting Bowtie2 alignment..."
bowtie2 -p 32 \
    -x ${GENOME_INDEX} \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    -I 0 -X 1000 \
    -1 analysis/trimmed/${sample}_R1_paired.fastq.gz \
    -2 analysis/trimmed/${sample}_R2_paired.fastq.gz \
    -S ${tmp_dir}/${sample}.sam \
    2> logs/aligned/${sample}.align.log

check_output ${tmp_dir}/${sample}.sam
log_message "Alignment completed successfully"

# Collect pre-filtering alignment statistics
log_message "Collecting pre-filtering statistics..."
samtools flagstat -@ 32 ${tmp_dir}/${sample}.sam > analysis/qc/${sample}_prefilter_flagstat.txt
log_message "Pre-filtering statistics saved to analysis/qc/${sample}_prefilter_flagstat.txt"

# Analyze mapping quality distribution
log_message "Checking quality score distribution of mapped reads..."
samtools view -@ 32 -q 1 ${tmp_dir}/${sample}.sam | \
    cut -f5 | sort -n | uniq -c > analysis/qc/${sample}_mapq_distribution.txt
log_message "MAPQ distribution saved to analysis/qc/${sample}_mapq_distribution.txt"

# Filter alignments with detailed logging
log_message "Converting SAM to BAM and filtering..."
log_message "Filtering criteria:"
log_message "  - Removing unmapped reads (-F 4)"
log_message "  - Minimum MAPQ score: 20 (-q 20)"

# Track read counts through filtering steps
total_reads=$(samtools view -c ${tmp_dir}/${sample}.sam)
log_message "Total reads before filtering: ${total_reads}"

# Filter and count reads in steps
samtools view -@ 32 -bh -F 4 ${tmp_dir}/${sample}.sam | \
    tee >(samtools view -c > ${tmp_dir}/mapped_count.txt) | \
    samtools view -@ 32 -bh -q 20 - > ${tmp_dir}/${sample}.bam

mapped_reads=$(cat ${tmp_dir}/mapped_count.txt)
final_reads=$(samtools view -c ${tmp_dir}/${sample}.bam)

# Log filtering statistics
log_message "Reads after removing unmapped: ${mapped_reads}"
log_message "Reads after MAPQ filtering: ${final_reads}"
log_message "Filtering summary:"
log_message "  - Lost in mapping: $((total_reads - mapped_reads)) reads ($(echo "scale=2; 100 * (1 - ${mapped_reads}/${total_reads})" | bc)%)"
log_message "  - Lost in MAPQ filter: $((mapped_reads - final_reads)) reads ($(echo "scale=2; 100 * (1 - ${final_reads}/${mapped_reads})" | bc)%)"

check_output ${tmp_dir}/${sample}.bam

# Collect pre-deduplication statistics
log_message "Collecting pre-duplicate marking statistics..."
samtools flagstat -@ 32 analysis/aligned/${sample}.sorted.bam > \
    analysis/qc/${sample}_prededup_flagstat.txt

# Mark duplicates and analyze results
log_message "Analyzing duplicate marking results..."
marked_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
dup_reads=$((final_reads - marked_reads))
dup_percent=$(echo "scale=2; 100 * ${dup_reads}/${final_reads}" | bc)

# Log duplication statistics
log_message "Duplication summary:"
log_message "  - Input reads: ${final_reads}"
log_message "  - Reads after duplicate marking: ${marked_reads}"
log_message "  - Duplicate reads: ${dup_reads} (${dup_percent}%)"

# Generate comprehensive alignment summary report
cat << EOF > analysis/qc/${sample}_alignment_summary.txt
=== Pipeline Summary for ${sample} ===
Initial reads from trimming: ${initial_reads}
Mapped reads: ${mapped_reads}
Reads passing MAPQ filter: ${final_reads}
Final reads after duplicate marking: ${marked_reads}

=== Loss at each step ===
Mapping loss: $((initial_reads - mapped_reads)) reads ($(echo "scale=2; 100 * (1 - ${mapped_reads}/${initial_reads})" | bc)%)
MAPQ filtering loss: $((mapped_reads - final_reads)) reads ($(echo "scale=2; 100 * (1 - ${final_reads}/${mapped_reads})" | bc)%)
Duplicate reads: ${dup_reads} (${dup_percent}%)

=== Final statistics ===
Overall retention rate: $(echo "scale=2; 100 * ${marked_reads}/${initial_reads}" | bc)%
EOF

log_message "Detailed alignment summary saved to analysis/qc/${sample}_alignment_summary.txt"

# Sort BAM file
log_message "Sorting BAM..."
samtools sort -@ 32 -m 2G ${tmp_dir}/${sample}.bam -o analysis/aligned/${sample}.sorted.bam
check_output analysis/aligned/${sample}.sorted.bam

# Clean up temporary files
log_message "Cleaning up temporary files..."
rm -rf ${tmp_dir}

# Calculate mapping statistics for raw alignments
get_mapping_stats analysis/aligned/${sample}.sorted.bam "raw"

# Index sorted BAM
log_message "Indexing BAM file..."
samtools index -@ 32 analysis/aligned/${sample}.sorted.bam

# Mark duplicates with Picard
# Uses more memory and removes optical duplicates
log_message "Marking duplicates..."
picard -Xmx48g MarkDuplicates \
    I=analysis/aligned/${sample}.sorted.bam \
    O=analysis/aligned/${sample}.dedup.bam \
    M=analysis/qc/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    VALIDATION_STRINGENCY=LENIENT

check_output analysis/aligned/${sample}.dedup.bam
check_output analysis/qc/${sample}_dup_metrics.txt

# Calculate mapping statistics for deduplicated alignments
get_mapping_stats analysis/aligned/${sample}.dedup.bam "dedup"

# Index final BAM
samtools index -@ 32 analysis/aligned/${sample}.dedup.bam

# Generate final alignment statistics
echo "Generating alignment statistics..."
samtools flagstat -@ 32 analysis/aligned/${sample}.dedup.bam > \
    analysis/qc/${sample}_flagstat.txt

# Generate fragment size distribution metrics and plot
echo "Generating fragment size distribution..."
picard -Xmx48g CollectInsertSizeMetrics \
    I=analysis/aligned/${sample}.dedup.bam \
    O=analysis/qc/${sample}_insert_size_metrics.txt \
    H=analysis/qc/${sample}_insert_size_histogram.pdf \
    M=0.5 \
    W=1000

# Clean up intermediate files to save space
if [[ -s analysis/aligned/${sample}.dedup.bam ]]; then
    rm -f analysis/aligned/${sample}.bam
    rm -f analysis/aligned/${sample}.sorted.bam
    rm -f analysis/aligned/${sample}.sorted.bam.bai
fi

log_message "Alignment pipeline completed successfully for ${sample}"