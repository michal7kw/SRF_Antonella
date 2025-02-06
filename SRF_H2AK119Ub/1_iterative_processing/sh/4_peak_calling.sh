#!/bin/bash

# This script performs peak calling and filtering on CUT&Tag data using MACS2
# It processes multiple samples in parallel using a SLURM array job
# For each sample, it:
# 1. Calls broad peaks using MACS2 with optimized parameters for CUT&Tag data
# 2. Filters peaks based on quality metrics and size
# 3. Removes blacklisted regions
# 4. Generates visualization-friendly peak files
# 5. Calculates and outputs comprehensive QC metrics

# Input files (per sample):
# - analysis/aligned/{sample}.dedup.bam - Deduplicated BAM file
# - analysis/aligned/{sample}.dedup.bam.bai - BAM index
# - hg38-blacklist.v2.bed - ENCODE blacklist regions

# Output files (per sample):
# - analysis/peaks/{sample}_broad_peaks.broadPeak - Raw MACS2 peaks
# - analysis/peaks/{sample}_broad_peaks_filtered.broadPeak - Quality filtered peaks
# - analysis/peaks/{sample}_broad_peaks_final.broadPeak - Blacklist filtered peaks
# - analysis/peaks/{sample}_broad_peaks_viz.broadPeak - Visualization-optimized peaks
# - analysis/qc/{sample}_metrics.csv - Peak calling QC metrics
# - analysis/peak_stats/{sample}_peak_stats.txt - Peak statistics

#SBATCH --job-name=4_peak_calling
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling.err"
#SBATCH --output="logs/4_peak_calling.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create output directories
mkdir -p analysis/{peaks,qc,peak_stats}

# Define sample names
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Verify input files exist
if [[ ! -f analysis/aligned/${sample}.dedup.bam ]] || [[ ! -f analysis/aligned/${sample}.dedup.bam.bai ]]; then
    echo "Error: Input BAM or index not found for ${sample}"
    exit 1
fi

echo "Processing ${sample} for broad peak calling..."

# Step 1: Initial MACS2 peak calling with parameters optimized for CUT&Tag
# --broad: Call broad peaks appropriate for histone modifications
# -q 0.05: Use FDR threshold of 5%
# --min-length 500: Minimum peak size of 500bp
# --buffer-size 10000: Increase buffer for large datasets
macs2 callpeak \
    -t analysis/aligned/${sample}.dedup.bam \
    -f BAMPE \
    -g hs \
    --broad \
    -n ${sample}_broad \
    --outdir analysis/peaks \
    -q 0.05 \
    --broad-cutoff 0.05 \
    --keep-dup 1 \
    --min-length 500 \
    --bdg \
    --buffer-size 10000 \
    --verbose 3

# Verify MACS2 output
if [ ! -s analysis/peaks/${sample}_broad_peaks.broadPeak ]; then
    echo "Error: MACS2 failed to produce output for ${sample}"
    exit 1
fi

# Step 2: Define file paths for processing stages
raw_peaks="analysis/peaks/${sample}_broad_peaks.broadPeak"
filtered_peaks="analysis/peaks/${sample}_broad_peaks_filtered.broadPeak"
final_peaks="analysis/peaks/${sample}_broad_peaks_final.broadPeak"
viz_peaks="analysis/peaks/${sample}_broad_peaks_viz.broadPeak"

# Step 3: Initial filtering and standardization of peak format
# Filters peaks based on:
# - Fold enrichment >= 3
# - Length between 500bp and 10kb
# - Standardizes peak names and scores for IGV compatibility
awk -v sample="$sample" 'BEGIN{OFS="\t"} {
    len = $3 - $2;
    # Store original values
    orig_score = $5;
    orig_fold = $7;
    orig_pval = $8;
    orig_qval = $9;
    
    # Basic filtering criteria
    if(orig_fold >= 3 && len >= 500 && len <= 10000) {
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
}' $raw_peaks > $filtered_peaks

# Add validation step to ensure BED format compliance
awk -v sample="$sample" '
BEGIN {
    valid=1;
    OFS="\t";
}
{
    # Validate each field
    if($1 !~ /^[a-zA-Z0-9_]+$/) {
        print "Invalid chromosome name at line " NR > "/dev/stderr";
        valid=0;
    }
    if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {
        print "Invalid coordinates at line " NR > "/dev/stderr";
        valid=0;
    }
    if($4 !~ /^[a-zA-Z0-9_]+$/) {
        print "Invalid peak name at line " NR > "/dev/stderr";
        valid=0;
    }
    if($5 !~ /^[0-9]+$/ || $5 < 0 || $5 > 1000) {
        print "Invalid score at line " NR > "/dev/stderr";
        valid=0;
    }
    if($6 !~ /^[.+-]$/) {
        print "Invalid strand at line " NR > "/dev/stderr";
        valid=0;
    }
    # If valid, print the line
    if(valid) print $0;
    valid=1;
}' $filtered_peaks > ${filtered_peaks}.tmp && \
mv ${filtered_peaks}.tmp $filtered_peaks

# Step 4: Remove peaks overlapping with blacklisted regions
bedtools intersect -v \
    -a $filtered_peaks \
    -b hg38-blacklist.v2.bed \
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
total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
reads_in_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
    -b $final_peaks -u | samtools view -c)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
peak_count=$(wc -l < $final_peaks)
mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' $final_peaks)
mean_fold=$(awk '{sum += $7} END {print sum/NR}' $final_peaks)

# Save detailed QC metrics to CSV
cat << EOF > analysis/qc/${sample}_metrics.csv
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
}' $final_peaks > analysis/peak_stats/${sample}_peak_stats.txt

echo "Completed processing ${sample}"
echo "Final peak count: $peak_count"
echo "FRiP score: $frip" 