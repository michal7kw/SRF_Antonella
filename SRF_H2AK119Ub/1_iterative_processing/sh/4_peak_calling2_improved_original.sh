#!/bin/bash

#SBATCH --job-name=4_peak_calling2_improved
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling2_improved.err"
#SBATCH --output="logs/4_peak_calling2_improved.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/{peaks2_improved,qc2_improved,peak_stats_improved}

# Define samples
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Check input files
if [[ ! -f analysis/aligned/${sample}.dedup.bam ]] || [[ ! -f analysis/aligned/${sample}.dedup.bam.bai ]]; then
    echo "Error: Input BAM or index not found for ${sample}"
    exit 1
fi

echo "Processing ${sample} for broad peak calling..."

# Step 1: Initial MACS2 peak calling
macs2 callpeak \
    -t analysis/aligned/${sample}.dedup.bam \
    -f BAMPE \
    -g hs \
    --broad \
    -n ${sample}_broad \
    --outdir analysis/peaks2_improved \
    -q 0.05 \
    --broad-cutoff 0.05 \
    --keep-dup auto \
    --min-length 500 \
    --bdg

# Check MACS2 output
if [ ! -s analysis/peaks2_improved/${sample}_broad_peaks.broadPeak ]; then
    echo "Error: MACS2 failed to produce output for ${sample}"
    exit 1
fi

# Step 2: Create temporary files for processing stages
raw_peaks="analysis/peaks2_improved/${sample}_broad_peaks.broadPeak"
filtered_peaks="analysis/peaks2_improved/${sample}_broad_peaks_filtered.broadPeak"
final_peaks="analysis/peaks2_improved/${sample}_broad_peaks_final.broadPeak"
viz_peaks="analysis/peaks2_improved/${sample}_broad_peaks_viz.broadPeak"

# Step 3: Initial filtering while preserving original scores
awk 'BEGIN{OFS="\t"} {
    len = $3 - $2;
    # Store original values
    orig_score = $5;
    orig_fold = $7;
    orig_pval = $8;
    orig_qval = $9;
    
    # Basic filtering criteria
    if(orig_fold >= 3 && len >= 500 && len <= 10000) {
        name = sprintf("%s_broad_peak_%d", "'${sample}'", NR);
        strand = ($6 == "+" || $6 == "-") ? $6 : ".";
        print $1, $2, $3, name, orig_score, strand, orig_fold, orig_pval, orig_qval;
    }
}' $raw_peaks > $filtered_peaks

# Step 4: Blacklist filtering
bedtools intersect -v \
    -a $filtered_peaks \
    -b hg38-blacklist.v2.bed \
    > $final_peaks

# Step 5: Create visualization-friendly version with normalized scores
awk 'BEGIN{OFS="\t"} {
    # Create normalized score for visualization (0-1000 range)
    viz_score = int($7 * 100);
    if(viz_score > 1000) viz_score = 1000;
    
    # Output same data but with normalized score for visualization
    print $1, $2, $3, $4, viz_score, $6, $7, $8, $9;
}' $final_peaks > $viz_peaks

# Step 6: Generate QC metrics
total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
reads_in_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
    -b $final_peaks -u | samtools view -c)
frip=$(echo "scale=4; $reads_in_peaks / $total_reads" | bc)
peak_count=$(wc -l < $final_peaks)
mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' $final_peaks)
mean_fold=$(awk '{sum += $7} END {print sum/NR}' $final_peaks)

# Save detailed QC metrics
cat << EOF > analysis/qc2_improved/${sample}_metrics.csv
Metric,Value
Sample,${sample}
Total_Reads,${total_reads}
Reads_in_Peaks,${reads_in_peaks}
FRiP,${frip}
Peak_Count,${peak_count}
Mean_Length,${mean_length}
Mean_Fold_Enrichment,${mean_fold}
EOF

# Step 7: Generate peak statistics
awk -v sample=$sample '
BEGIN {OFS="\t"}
{
    len = $3 - $2;
    fold = $7;
    print len, fold;
}' $final_peaks > analysis/peak_stats_improved/${sample}_peak_stats.txt

echo "Completed processing ${sample}"
echo "Final peak count: $peak_count"
echo "FRiP score: $frip" 