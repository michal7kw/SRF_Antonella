#!/bin/bash

# Usage: bash 5_peak_calling_strict.sh <sample_index>
# Example: bash 5_peak_calling_strict.sh 0  (for GFP_1)
# Note: This script processes one sample at a time. Run multiple instances for parallel processing.

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
# It processes multiple samples sequentially based on command-line input.
# For each sample, it:
# 1. Calls broad peaks using MACS2 with optimized parameters for CUT&Tag data
# 2. Filters peaks based on quality metrics and size
# 3. Removes blacklisted regions
# 4. Generates visualization-friendly peak files
# 5. Calculates and outputs comprehensive QC metrics

# Input files (per sample):
# - analysis/3_alignment/{sample}.dedup.bam - Deduplicated BAM file
# - analysis/3_alignment/{sample}.dedup.bam.bai - BAM index
# - ../COMMON_DATA/hg38-blacklist.v2.bed - ENCODE blacklist regions

# Output files (per sample):
# - analysis/5_peak_calling_strict/{sample}_broad_peaks.broadPeak - Raw MACS2 peaks
# - analysis/5_peak_calling_strict/{sample}_broad_peaks_filtered.broadPeak - Quality filtered peaks
# - analysis/5_peak_calling_strict/{sample}_broad_peaks_final.broadPeak - Blacklist filtered peaks
# - analysis/5_peak_calling_strict/{sample}_broad_peaks_viz.broadPeak - Visualization-optimized peaks
# - analysis/5_peak_calling_strict/{sample}_metrics.csv - Peak calling QC metrics
# - analysis/5_peak_calling_strict/{sample}_peak_stats.txt - Peak statistics

set -e  # Exit on error
set -u  # Exit on undefined variable

# Check for command-line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <sample_index>"
    echo "Provide the index (0-5) of the sample to process."
    exit 1
fi
SAMPLE_INDEX=$1

# Activate conda environment with required tools
# Ensure the 'snakemake' conda environment is active before running this script
# Example: conda activate snakemake
# source /opt/common/tools/ric.cosr/miniconda3/bin/activate # Removed cluster-specific path
# conda activate snakemake # Assuming environment is already active

# Create a unique temporary directory for this job run
# Using a timestamp and the sample index for uniqueness locally
TIMESTAMP=$(date +%Y%m%d%H%M%S)
TMPDIR="tmp_strict/${SAMPLE_INDEX}_${TIMESTAMP}"
mkdir -p "$TMPDIR"
# chmod 777 "$TMPDIR" # Permissions might not be needed locally or could be adjusted
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
# unset PYTHONPATH # Keep if needed locally, MACS2 export above should override
unset PYTHONHOME
unset PYTHONEXECUTABLE

# Define working directory (assuming script is run from 1_iterative_processing)
WORKDIR="." # Current directory
cd $WORKDIR

# Define output directory (make it specific to strict run)
OUTPUT_DIR="analysis/5_peak_calling_strict"

# Create output directories
mkdir -p ${OUTPUT_DIR}
mkdir -p logs # Create logs directory if it doesn't exist

# Define sample names
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# Validate sample index
if [[ $SAMPLE_INDEX -lt 0 ]] || [[ $SAMPLE_INDEX -ge ${#samples[@]} ]]; then
    echo "Error: Invalid sample index $SAMPLE_INDEX. Must be between 0 and $((${#samples[@]} - 1))."
    exit 1
fi
sample=${samples[$SAMPLE_INDEX]}

# Verify input files exist
INPUT_BAM="analysis/3_alignment/${sample}.dedup.bam"
INPUT_BAI="analysis/3_alignment/${sample}.dedup.bam.bai"
BLACKLIST_BED="../COMMON_DATA/hg38-blacklist.v2.bed" # Relative path

if [[ ! -f $INPUT_BAM ]]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi
if [[ ! -f $INPUT_BAI ]]; then
    echo "Error: Input BAM index not found: $INPUT_BAI"
    exit 1
fi
if [[ ! -f $BLACKLIST_BED ]]; then
    echo "Error: Blacklist BED not found: $BLACKLIST_BED"
    exit 1
fi


echo "Processing ${sample} (Index: ${SAMPLE_INDEX}) for *strict* broad peak calling..."

# Create sample-specific temp directory
SAMPLE_TMP="${TMPDIR}/${sample}"
mkdir -p "$SAMPLE_TMP"
# chmod 777 "$SAMPLE_TMP" # Permissions might not be needed locally

# Step 0: Standardize chromosome names in BAM file
echo "Standardizing chromosome names for ${sample}..."
samtools view -H $INPUT_BAM | \
    sed 's/SN:\([0-9XY]\)/SN:chr\1/' | \
    sed 's/SN:MT/SN:chrM/' > "${SAMPLE_TMP}/${sample}_header.sam"

samtools reheader "${SAMPLE_TMP}/${sample}_header.sam" $INPUT_BAM > "${SAMPLE_TMP}/${sample}.chr.bam"
samtools index "${SAMPLE_TMP}/${sample}.chr.bam"

# Step 1: Initial MACS2 peak calling with parameters optimized for CUT&Tag (Strict)
echo "Running MACS2 peak calling for ${sample}..."
# Create MACS2 temp directory with unique path
MACS2_TMP="${SAMPLE_TMP}/macs2_tmp_${RANDOM}"
mkdir -p "$MACS2_TMP"
# chmod 777 "$MACS2_TMP" # Permissions might not be needed locally

# Run MACS2 with optimized settings (Removed env -i)
# Ensure necessary variables are exported (done above)
macs2 callpeak \
    -t "${SAMPLE_TMP}/${sample}.chr.bam" \
    -f BAMPE \
    -g hs \
    --broad \
    -n ${sample}_broad \
    --outdir ${OUTPUT_DIR} \
    -q 0.01 \
    --broad-cutoff 0.01 \
    --keep-dup 1 \
    --min-length 1000 \
    --bdg \
    --buffer-size 1000000 \
    --verbose 3 \
    --tempdir "$MACS2_TMP" \
    --nomodel \
    --nolambda \
    --slocal 1000 \
    --llocal 10000

# Verify MACS2 output
if [ ! -s ${OUTPUT_DIR}/${sample}_broad_peaks.broadPeak ]; then
    echo "Error: MACS2 failed to produce output for ${sample}"
    # Consider removing partial temp dir on error
    # rm -rf "$MACS2_TMP"
    # rm -rf "$SAMPLE_TMP"
    # rm -rf "$TMPDIR" # Careful with this if running in parallel
    exit 1
fi

# Step 2: Define file paths for processing stages
raw_peaks="${OUTPUT_DIR}/${sample}_broad_peaks.broadPeak"
filtered_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_filtered.broadPeak"
final_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_final.broadPeak"
viz_peaks="${OUTPUT_DIR}/${sample}_broad_peaks_viz.broadPeak"
metrics_csv="${OUTPUT_DIR}/${sample}_metrics.csv"
stats_txt="${OUTPUT_DIR}/${sample}_peak_stats.txt"

# Step 3: Initial filtering and standardization of peak format
echo "Filtering and standardizing peaks for ${sample}..."
# Fix chromosome names first
awk 'BEGIN{OFS="\t"} {
    # Fix chromosome names if needed
    if ($1 !~ /^chr/) {
        if ($1 ~ /^[0-9XY]+$/) {
            $1 = "chr" $1
        } else if ($1 == "MT") {
            $1 = "chrM"
        } else if ($1 ~ /^GL/) {
            # Skip scaffolds/patches or map them if necessary
            # For simplicity, skipping them here
            next
        } else if ($1 ~ /^KI/) {
             # Skip scaffolds/patches
             next
        } else {
            # Handle other unexpected chromosome names if necessary
            print "Warning: Skipping unexpected chromosome format: " $1 > "/dev/stderr"
            next
        }
    }
    print $0
}' $raw_peaks > "${SAMPLE_TMP}/fixed_chrom.bed"

# Apply filtering criteria (Strict: fold>=5, len 1000-50000)
awk -v sample="$sample" 'BEGIN{OFS="\t"} {
    len = $3 - $2;
    # Store original values
    orig_score = $5;
    orig_fold = $7;
    orig_pval = $8;
    orig_qval = $9;

    # Strict filtering criteria
    if(orig_fold >= 5 && len >= 1000 && len <= 50000) {
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

# Validate BED format (simplified)
echo "Validating filtered BED format..."
awk 'BEGIN{OFS="\t"} {
    if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 < 0 || $3 < $2) {
        print "Invalid coordinates at line " NR ": " $0 > "/dev/stderr"; next
    }
    if($5 !~ /^[0-9]+$/ || $5 < 0 || $5 > 1000) {
        print "Invalid score at line " NR ": " $0 > "/dev/stderr"; next
    }
    if($6 != "." && $6 != "+" && $6 != "-") {
         print "Invalid strand at line " NR ": " $0 > "/dev/stderr"; next
    }
    print $0
}' $filtered_peaks > ${filtered_peaks}.tmp && \
mv ${filtered_peaks}.tmp $filtered_peaks || { echo "BED validation failed"; exit 1; }

# Step 4: Remove peaks overlapping with blacklisted regions
echo "Removing blacklisted regions..."
bedtools intersect -v \
    -a $filtered_peaks \
    -b $BLACKLIST_BED \
    -nonamecheck \
    > $final_peaks

# Step 5: Create IGV-friendly peak file with normalized scores
echo "Creating visualization file..."
awk 'BEGIN{OFS="\t"} {
    # Create normalized score for visualization (0-1000 range based on fold enrichment)
    # Cap fold enrichment at a reasonable max (e.g., 50) for scaling
    max_fold = 50;
    scaled_fold = ($7 > max_fold) ? max_fold : $7;
    viz_score = int((scaled_fold / max_fold) * 1000);
    if(viz_score > 1000) viz_score = 1000;
    if(viz_score < 0) viz_score = 0; # Ensure non-negative

    # Output standard BED6 + 3 additional columns (signalValue, pValue, qValue)
    # Using normalized score as the main score (column 5)
    print $1, $2, $3, $4, viz_score, $6, $7, $8, $9;
}' $final_peaks > $viz_peaks


# Step 6: Calculate and output QC metrics
echo "Calculating QC metrics..."
# Use original non-chr-prefix BAM for total reads count
total_reads=$(samtools view -c $INPUT_BAM)

# Use the temporary chr-prefix BAM for FRiP calculation as peaks have chr prefix
TEMP_CHR_BAM="${SAMPLE_TMP}/${sample}.chr.bam"
if [ -f "$TEMP_CHR_BAM" ]; then
    # Use bedtools coverage for more accurate FRiP
    # Count reads where at least 1bp overlaps a peak
    reads_in_peaks=$(bedtools intersect -u -a "$TEMP_CHR_BAM" -b $final_peaks -f 0.5 -F 0.5 -e | samtools view -c -)

    if [ "$total_reads" -gt 0 ]; then
        frip=$(awk -v r=$reads_in_peaks -v t=$total_reads 'BEGIN {printf "%.6f", r/t}')
    else
        frip=0
    fi
else
    echo "Warning: ${TEMP_CHR_BAM} not found for FRiP calculation. Skipping FRiP."
    reads_in_peaks="NA"
    frip="NA"
fi


peak_count=$(wc -l < $final_peaks | awk '{print $1}') # Ensure only number is captured
if [ "$peak_count" -gt 0 ]; then
    mean_length=$(awk '{sum += $3-$2} END {printf "%.2f", sum/NR}' $final_peaks)
    mean_fold=$(awk '{sum += $7} END {printf "%.2f", sum/NR}' $final_peaks)
else
    mean_length=0
    mean_fold=0
fi

# Save detailed QC metrics to CSV
cat << EOF > $metrics_csv
Metric,Value
Sample,${sample}
Total_Reads,${total_reads}
Reads_in_Peaks,${reads_in_peaks}
FRiP,${frip}
Final_Peak_Count,${peak_count}
Mean_Peak_Length,${mean_length}
Mean_Fold_Enrichment,${mean_fold}
EOF

# Step 7: Generate detailed peak statistics for downstream analysis
echo "Generating peak statistics file..."
if [ "$peak_count" -gt 0 ]; then
    awk 'BEGIN {OFS="\t"; print "Length", "FoldEnrichment"} { print $3-$2, $7 }' $final_peaks > $stats_txt
else
    echo -e "Length\tFoldEnrichment" > $stats_txt # Create empty file with header
fi

echo "Completed *strict* processing for ${sample}"
echo "Final peak count: $peak_count"
echo "FRiP score: $frip"
echo "Output files are in ${OUTPUT_DIR}"
echo "Temporary files are in ${SAMPLE_TMP}"


# Clean up temporary files at the very end
# echo "Cleaning up temporary directory: ${SAMPLE_TMP}"
# rm -rf "${SAMPLE_TMP}"
# Consider leaving tmp dir for debugging or add option to keep/remove
echo "Note: Temporary directory ${SAMPLE_TMP} was not removed."

echo "Script finished successfully for ${sample}."