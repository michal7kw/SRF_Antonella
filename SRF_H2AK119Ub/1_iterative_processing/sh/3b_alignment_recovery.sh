#!/bin/bash
#SBATCH --job-name=3b_alignment_recovery
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/3b_alignment_recovery_%a.err"
#SBATCH --output="logs/3b_alignment_recovery_%a.out"

# This script continues the alignment pipeline from intermediate results
# It checks for existing files and only runs the remaining necessary steps

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

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
mkdir -p analysis/{aligned,qc} logs/aligned

# Define samples array and get current sample
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

log_message "Processing sample: ${sample}"

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
    total_reads=$(samtools view -c $bam)
    mapped_reads=$(samtools view -c -F 4 $bam)
    mapping_rate=$(echo "scale=4; $mapped_reads / $total_reads * 100" | bc)
    
    echo "${prefix},${total_reads},${mapped_reads},${mapping_rate}" > analysis/qc/${sample}_mapping_stats.csv
}

# Check for existing intermediate BAM
tmp_dir="analysis/aligned/tmp_${sample}"
if [[ ! -s "${tmp_dir}/${sample}.bam" ]]; then
    log_message "ERROR: Intermediate BAM file not found at ${tmp_dir}/${sample}.bam"
    exit 1
fi

# Continue from intermediate BAM file
log_message "Continuing from existing intermediate BAM file..."

# Get initial read counts from trimming metrics
if [[ -f analysis/qc/${sample}_trimming_metrics.json ]]; then
    trim_metrics=$(cat analysis/qc/${sample}_trimming_metrics.json)
    initial_reads=$(echo $trim_metrics | jq -r '.clean_r1_reads')
    log_message "Found trimming metrics: ${initial_reads} input reads"
fi

# Sort and filter BAM
log_message "Sorting and filtering BAM file..."
samtools view -@ 32 -bh -q 20 "${tmp_dir}/${sample}.bam" | \
    samtools sort -@ 32 -m 2G - -o analysis/aligned/${sample}.sorted.bam

check_output analysis/aligned/${sample}.sorted.bam

# Always add (or replace) Read Group (RG) tags to ensure complete information
log_message "Adding/replace Read Group information in ${sample}.sorted.bam..."
picard -Xmx16g AddOrReplaceReadGroups \
    I=analysis/aligned/${sample}.sorted.bam \
    O=analysis/aligned/${sample}.sorted.rg.bam \
    RGID=RG_${sample} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${sample} \
    VALIDATION_STRINGENCY=LENIENT

check_output analysis/aligned/${sample}.sorted.rg.bam
sorted_input="analysis/aligned/${sample}.sorted.rg.bam"

# Get read counts for reporting
final_reads=$(samtools view -c "$sorted_input")
log_message "Reads after MAPQ filtering: ${final_reads}"

# Collect pre-deduplication statistics
log_message "Collecting pre-duplicate marking statistics..."
samtools flagstat -@ 32 "$sorted_input" > \
    analysis/qc/${sample}_prededup_flagstat.txt

# Mark duplicates
log_message "Marking duplicates..."
picard -Xmx48g MarkDuplicates \
    I="$sorted_input" \
    O=analysis/aligned/${sample}.dedup.bam \
    M=analysis/qc/${sample}_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    VALIDATION_STRINGENCY=LENIENT

check_output analysis/aligned/${sample}.dedup.bam
check_output analysis/qc/${sample}_dup_metrics.txt

# Calculate final statistics
marked_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
dup_reads=$((final_reads - marked_reads))
dup_percent=$(echo "scale=2; 100 * ${dup_reads}/${final_reads}" | bc)

# Generate comprehensive alignment summary
cat << EOF > analysis/qc/${sample}_alignment_summary.txt
=== Pipeline Summary for ${sample} ===
Initial reads from trimming: ${initial_reads}
Final reads before deduplication: ${final_reads}
Final reads after duplicate marking: ${marked_reads}

=== Duplication Summary ===
Duplicate reads: ${dup_reads} (${dup_percent}%)

=== Final statistics ===
Overall retention rate: $(echo "scale=2; 100 * ${marked_reads}/${initial_reads}" | bc)%
EOF

# Index final BAM
samtools index -@ 32 analysis/aligned/${sample}.dedup.bam

# Generate final alignment statistics
log_message "Generating final alignment statistics..."
samtools flagstat -@ 32 analysis/aligned/${sample}.dedup.bam > \
    analysis/qc/${sample}_flagstat.txt

# Generate fragment size distribution
log_message "Generating fragment size distribution..."
picard -Xmx48g CollectInsertSizeMetrics \
    I=analysis/aligned/${sample}.dedup.bam \
    O=analysis/qc/${sample}_insert_size_metrics.txt \
    H=analysis/qc/${sample}_insert_size_histogram.pdf \
    M=0.5 \
    W=1000

# Clean up intermediate files if final BAM exists
if [[ -s analysis/aligned/${sample}.dedup.bam ]]; then
    log_message "Cleaning up intermediate files..."
    rm -f analysis/aligned/${sample}.sorted.bam
    rm -f analysis/aligned/${sample}.sorted.rg.bam
    rm -f analysis/aligned/${sample}.sorted.bam.bai
    rm -rf "${tmp_dir}"
fi

log_message "Recovery pipeline completed successfully for ${sample}" 