#!/bin/bash
#SBATCH --job-name=3_alignment
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3_alignment/3_alignment_%a.err"
#SBATCH --output="logs/3_alignment/3_alignment_%a.out"

# Exit on error, but we want to see the error message
set -e
set -u
set -o pipefail

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check if a command exists
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

# Print conda environment info
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

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/{aligned,qc} logs/aligned || { log_message "ERROR: Failed to create directories"; exit 1; }

# Function to check if output exists and is not empty
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

# Define samples
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

# Create temporary directory for this sample
tmp_dir="analysis/aligned/tmp_${sample}"
log_message "Creating temporary directory: ${tmp_dir}"
mkdir -p ${tmp_dir} || { log_message "ERROR: Failed to create temporary directory"; exit 1; }

# Run alignment and save as SAM
log_message "Starting Bowtie2 alignment..."
bowtie2 -p 32 \
    -x ${GENOME_INDEX} \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    --no-overlap --no-dovetail \
    -I 10 -X 700 \
    -1 analysis/trimmed/${sample}_R1_paired.fastq.gz \
    -2 analysis/trimmed/${sample}_R2_paired.fastq.gz \
    -S ${tmp_dir}/${sample}.sam \
    2> logs/aligned/${sample}.align.log

check_output ${tmp_dir}/${sample}.sam
log_message "Alignment completed successfully"

# Convert to BAM, filter, and sort in separate steps
log_message "Converting SAM to BAM and filtering..."
samtools view -@ 32 -bh -F 4 -q 30 ${tmp_dir}/${sample}.sam > ${tmp_dir}/${sample}.bam
check_output ${tmp_dir}/${sample}.bam

log_message "Sorting BAM..."
samtools sort -@ 32 -m 2G ${tmp_dir}/${sample}.bam -o analysis/aligned/${sample}.sorted.bam
check_output analysis/aligned/${sample}.sorted.bam

# Remove temporary files
log_message "Cleaning up temporary files..."
rm -rf ${tmp_dir}

get_mapping_stats analysis/aligned/${sample}.sorted.bam "raw"

# Index BAM
log_message "Indexing BAM file..."
samtools index -@ 32 analysis/aligned/${sample}.sorted.bam

# Mark duplicates with more memory and remove optical duplicates
log_message "Marking duplicates..."
picard -Xmx48g MarkDuplicates \
    I=analysis/aligned/${sample}.sorted.bam \
    O=analysis/aligned/${sample}.dedup.bam \
    M=analysis/qc/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_SEQUENCING_DUPLICATES=true

check_output analysis/aligned/${sample}.dedup.bam
check_output analysis/qc/${sample}_dup_metrics.txt

get_mapping_stats analysis/aligned/${sample}.dedup.bam "dedup"

# Index final BAM
samtools index -@ 32 analysis/aligned/${sample}.dedup.bam

# Generate comprehensive alignment statistics
echo "Generating alignment statistics..."
samtools flagstat -@ 32 analysis/aligned/${sample}.dedup.bam > \
    analysis/qc/${sample}_flagstat.txt

# Generate fragment size distribution
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