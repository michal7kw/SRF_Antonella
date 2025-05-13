#!/bin/bash

#SBATCH --job-name=5_peak_calling
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/5_peak_calling.err"
#SBATCH --output="logs/5_peak_calling.out"

# Documentation:
# This script performs peak calling and filtering on CUT&Tag data using MACS2
# It processes multiple samples in parallel using a SLURM array job
# For each sample, it:
# 1. Calls broad peaks using MACS2 with optimized parameters for CUT&Tag data
# 2. Filters peaks based on quality metrics and size
# 3. Removes blacklisted regions
# 4. Generates visualization-friendly peak files
# 5. Calculates and outputs comprehensive QC metrics

# Input files (per sample):
# - analysis/3_alignment/{sample}.dedup.bam - Deduplicated BAM file
# - analysis/3_alignment/{sample}.dedup.bam.bai - BAM index
# - hg38-blacklist.v2.bed - ENCODE blacklist regions

# Output files (per sample):
# - {sample}_broad_peaks.broadPeak - Raw MACS2 peaks
# - {sample}_broad_peaks_filtered.broadPeak - Quality filtered peaks
# - {sample}_broad_peaks_final.broadPeak - Blacklist filtered peaks
# - {sample}_broad_peaks_viz.broadPeak - Visualization-optimized peaks
# - {sample}_metrics.csv - Peak calling QC metrics
# - {sample}_peak_stats.txt - Peak statistics

set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create a unique temporary directory for this job
TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/tmp/${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
chmod 777 "$TMPDIR"
export TMPDIR

# Set Python environment variables for large data handling
export PYTHONPATH=""
export MACS2_LARGE_DATA=1
export PYTHONMALLOC=malloc
export PYTHONWARNINGS=ignore
export PYTHONHASHSEED=0
export PYTHONDONTWRITEBYTECODE=1

# Clean any existing Python environment variables
unset PYTHONSTARTUP
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONEXECUTABLE

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling"

# Create output directories
mkdir -p ${OUTPUT_DIR}

# Define sample names
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Verify input files exist
if [[ ! -f analysis/3_alignment/${sample}.dedup.bam ]] || [[ ! -f analysis/3_alignment/${sample}.dedup.bam.bai ]]; then
    echo "Error: Input BAM or index not found for ${sample}"
    exit 1
fi

echo "Processing ${sample} for broad peak calling..."

# Create sample-specific temp directory
SAMPLE_TMP="${TMPDIR}/${sample}"
mkdir -p "$SAMPLE_TMP"
chmod 777 "$SAMPLE_TMP"

# Step 0: Standardize chromosome names in BAM file
echo "Standardizing chromosome names for ${sample}..."
samtools view -H analysis/3_alignment/${sample}.dedup.bam | \
    sed 's/SN:\([0-9XY]\)/SN:chr\1/' | \
    sed 's/SN:MT/SN:chrM/' > "${SAMPLE_TMP}/${sample}_header.sam"

samtools reheader "${SAMPLE_TMP}/${sample}_header.sam" analysis/3_alignment/${sample}.dedup.bam > "${SAMPLE_TMP}/${sample}.chr.bam"
samtools index "${SAMPLE_TMP}/${sample}.chr.bam"

# Step 1: Initial MACS2 peak calling with parameters optimized for CUT&Tag
echo "Running MACS2 peak calling for ${sample}..."
# Create MACS2 temp directory with unique path
MACS2_TMP="${SAMPLE_TMP}/macs2_tmp_${RANDOM}"
mkdir -p "$MACS2_TMP"
chmod 777 "$MACS2_TMP"

# Run MACS2 with the specified parameters
env -i \
    HOME="$HOME" \
    PATH="$PATH" \
    TMPDIR="$MACS2_TMP" \
    PYTHONPATH="" \
    MACS2_LARGE_DATA=1 \
    PYTHONMALLOC=malloc \
    macs2 callpeak \
    -t "${SAMPLE_TMP}/${sample}.chr.bam" \
    -f BAMPE \
    --broad \
    -g hs \
    --broad-cutoff 0.05 \
    --nomodel \
    --extsize 200 \
    -q 0.05 \
    -n ${sample}_broad \
    --outdir ${OUTPUT_DIR} \
    --keep-dup all \
    --bdg \
    --buffer-size 1000000 \
    --verbose 3 \
    --tempdir "$MACS2_TMP"

# Verify MACS2 output
if [ ! -s ${OUTPUT_DIR}/${sample}_broad_peaks.broadPeak ]; then
    echo "Error: MACS2 failed to produce output for ${sample}"
    exit 1
fi

# Step 2: Define file paths for processing stages
raw_peaks="${OUTPUT_DIR}/${sample}_broad_peaks.broadPeak"
filtered_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_filtered.broadPeak"
final_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_final.broadPeak"
viz_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_viz.broadPeak"

# Step 3: Initial filtering and standardization of peak format
# Don't delete the temporary directory yet - we need it for FRiP calculation later
# rm -rf "${SAMPLE_TMP}"

# Fix chromosome names to ensure they match the expected format
# First, create a version with standardized chromosome names
awk 'BEGIN{OFS="\t"} {
    # Fix chromosome names if needed
    if ($1 !~ /^chr/) {
        if ($1 ~ /^[0-9XY]+$/) {
            $1 = "chr" $1
        } else if ($1 == "MT") {
            $1 = "chrM"
        } else if ($1 ~ /^GL/) {
            $1 = "chr" $1
        } else if ($1 ~ /^KI/) {
            $1 = "chr" $1
        }
    }
    print $0
}' $raw_peaks > "${SAMPLE_TMP}/fixed_chrom.bed"

# Now apply the filtering criteria with more lenient parameters for H2AK119Ub
awk -v sample="$sample" 'BEGIN{OFS="\t"} {
    len = $3 - $2;
    # Store original values
    orig_score = $5;
    orig_fold = $7;
    orig_pval = $8;
    orig_qval = $9;
    
    # More lenient filtering criteria for H2AK119Ub
    # Reduced fold enrichment threshold and broader size range
    if(orig_fold >= 2 && len >= 500 && len <= 100000) {
        # Ensure name follows strict BED format (no special characters)
        name = sprintf("%s_peak_%d", sample, NR);
        name = gensub(/[^a-zA-Z0-9_]/, "_", "g", name);  # Replace any invalid chars with underscore
        
        # Ensure score is integer between 0-1000
        score = int(orig_score);
        if(score > 1000) score = 1000;
        if(score < 0) score = 0;
        
        # Ensure strand is valid
        strand = ($6 == "+" || $6 == "-") ? $6 : ".";
        
        print $1, $2, $3, name, score, strand, orig_fold, orig_pval, orig_qval;
    }
}' "${SAMPLE_TMP}/fixed_chrom.bed" > $filtered_peaks

# Add validation step to ensure BED format compliance
# Simplified to avoid chromosome name validation errors
awk -v sample="$sample" '
BEGIN {
    OFS="\t";
}
{
    # Only validate numeric fields and output all lines
    if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {
        print "Invalid coordinates at line " NR > "/dev/stderr";
    }
    else if($5 !~ /^[0-9]+$/ || $5 < 0 || $5 > 1000) {
        print "Invalid score at line " NR > "/dev/stderr";
    }
    else {
        print $0;
    }
}' $filtered_peaks > ${filtered_peaks}.tmp && \
mv ${filtered_peaks}.tmp $filtered_peaks

# Step 4: Remove peaks overlapping with blacklisted regions
bedtools intersect -v \
    -a $filtered_peaks \
    -b /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/hg38-blacklist.v2.bed \
    -nonamecheck \
    > $final_peaks

# Step 5: Create IGV-friendly peak file with normalized scores
awk 'BEGIN{OFS="\t"} {
    # Create normalized score for visualization (0-1000 range)
    viz_score = int($7 * 100);
    if(viz_score > 1000) viz_score = 1000;
    
    # Output same data but with normalized score for visualization
    print $1, $2, $3, $4, viz_score, $6, $7, $8, $9;
}' $final_peaks > $viz_peaks

# Step 6: Calculate and output QC metrics
# - Total reads in sample
# - Reads in peaks
# - Fraction of reads in peaks (FRiP)
# - Peak count and characteristics
total_reads=$(samtools view -c analysis/3_alignment/${sample}.dedup.bam)

# Fix FRiP calculation to properly count reads in peaks
# Make sure we're using the chr.bam file that still exists
if [ -f "${SAMPLE_TMP}/${sample}.chr.bam" ]; then
    reads_in_peaks=$(bedtools intersect -a "${SAMPLE_TMP}/${sample}.chr.bam" \
        -b $final_peaks -c -bed | awk '$NF > 0' | wc -l)
    
    # Calculate FRiP score with proper precision
    frip=$(awk -v r=$reads_in_peaks -v t=$total_reads 'BEGIN {printf "%.6f", r/t}')
else
    # Fallback if the chr.bam file was deleted
    echo "Warning: ${SAMPLE_TMP}/${sample}.chr.bam not found for FRiP calculation. Using original BAM."
    reads_in_peaks=$(bedtools intersect -a analysis/3_alignment/${sample}.dedup.bam \
        -b $final_peaks -c -bed | awk '$NF > 0' | wc -l)
    
    # Calculate FRiP score with proper precision
    frip=$(awk -v r=$reads_in_peaks -v t=$total_reads 'BEGIN {printf "%.6f", r/t}')
fi

peak_count=$(wc -l < $final_peaks)
mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' $final_peaks)
mean_fold=$(awk '{sum += $7} END {print sum/NR}' $final_peaks)

# Save detailed QC metrics to CSV
cat << EOF > ${OUTPUT_DIR}/${sample}_metrics.csv
Metric,Value
Sample,${sample}
Total_Reads,${total_reads}
Reads_in_Peaks,${reads_in_peaks}
FRiP,${frip}
Peak_Count,${peak_count}
Mean_Length,${mean_length}
Mean_Fold_Enrichment,${mean_fold}
EOF

# Step 7: Generate detailed peak statistics for downstream analysis
awk -v sample=$sample '
BEGIN {OFS="\t"}
{
    len = $3 - $2;
    fold = $7;
    print len, fold;
}' $final_peaks > ${OUTPUT_DIR}/${sample}_peak_stats.txt

echo "Completed processing ${sample}"
echo "Final peak count: $peak_count"
echo "FRiP score: $frip"

# Clean up temporary files at the very end
rm -rf "${SAMPLE_TMP}"