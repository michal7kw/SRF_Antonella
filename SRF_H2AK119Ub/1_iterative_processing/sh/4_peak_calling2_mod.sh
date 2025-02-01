#!/bin/bash

#SBATCH --job-name=4_peak_calling2_mod
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling/4_peak_calling2_mod.err"
#SBATCH --output="logs/4_peak_calling/4_peak_calling2_mod.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/{peaks2_mod,qc2_mod} 
#logs/peaks2_mod

# Define samples
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
analysis_types=(broad)

# Since we only have one analysis type (broad), we don't need to calculate analysis_idx
sample_idx=$SLURM_ARRAY_TASK_ID
sample=${samples[$sample_idx]}
analysis_type=${analysis_types[0]}

# Function to check if input/output exists and is not empty
check_output() {
    local file=$1
    if [[ ! -s $file ]]; then
        echo "Error: Output file $file is empty or does not exist"
        exit 1
    fi
}

# Check if input BAM exists and is properly indexed
if [[ ! -f analysis/aligned/${sample}.dedup.bam ]] || [[ ! -f analysis/aligned/${sample}.dedup.bam.bai ]]; then
    echo "Error: Input BAM or index not found for ${sample}"
    exit 1
fi

echo "Processing ${sample} for ${analysis_type} peak calling..."

# Peak calling based on analysis type
if [[ $analysis_type == "broad" ]]; then
    echo "Calling broad peaks for ${sample} with stringent parameters..."

    # Modified MACS2 command to handle large datasets better
    macs2 callpeak \
        -t analysis/aligned/${sample}.dedup.bam \
        -f BAMPE \
        -g hs \
        --broad \
        -n ${sample}_broad \
        --outdir analysis/peaks2_mod \
        -q 0.05 \
        --broad-cutoff 0.05 \
        --keep-dup auto \
        --min-length 500 \
        --bdg

    # Add error checking for MACS2
    if [ $? -ne 0 ]; then
        echo "MACS2 failed for ${sample}"
        exit 1
    fi

    check_output analysis/peaks2_mod/${sample}_broad_peaks.broadPeak

    # Generate statistics about unfiltered peaks
    echo "=== Peak Statistics Before Filtering ==="
    echo "Total peaks: $(wc -l < analysis/peaks2_mod/${sample}_broad_peaks.broadPeak)"
    awk '{sum+=$7; if($7>max) max=$7} END {print "Average fold enrichment: "sum/NR"\nMax fold enrichment: "max}' analysis/peaks2_mod/${sample}_broad_peaks.broadPeak
    awk '{len=$3-$2; sum+=len} END {print "Average peak length: "sum/NR}' analysis/peaks2_mod/${sample}_broad_peaks.broadPeak

    # Post-filtering with more flexible criteria and ensuring IGV compatibility
    # .broadPeak format (BED6+3):
    # 1. chrom
    # 2. chromStart
    # 3. chromEnd
    # 4. name
    # 5. score (0-1000)
    # 6. strand (+/- or .)
    # 7. signalValue (fold enrichment)
    # 8. pValue (-log10)
    # 9. qValue (-log10)
    awk 'BEGIN{OFS="\t"} {
        name=sprintf("%s_broad_peak_%d", "'${sample}'", NR);
        len=$3-$2;
        # Normalize score to 0-1000 range for IGV
        score=int($7 * 100);
        if(score > 1000) score=1000;
        # Filter based on:
        # 1. Fold enrichment >= 3
        # 2. Peak length between 500bp and 10000bp (typical for broad marks)
        if($7 >= 3 && len >= 500 && len <= 10000) {
            # Ensure proper strand format
            strand=($6 == "+" || $6 == "-") ? $6 : ".";
            # Print in proper .broadPeak format
            print $1, $2, $3, name, score, strand, $7, $8, $9
        }
    }' analysis/peaks2_mod/${sample}_broad_peaks.broadPeak > analysis/peaks2_mod/${sample}_broad_peaks_filtered.broadPeak

    # Generate statistics about filtered peaks
    echo "=== Peak Statistics After Filtering ==="
    echo "Remaining peaks: $(wc -l < analysis/peaks2_mod/${sample}_broad_peaks_filtered.broadPeak)"
    awk '{sum+=$7} END {print "Average fold enrichment after filtering: "sum/NR}' analysis/peaks2_mod/${sample}_broad_peaks_filtered.broadPeak
    awk '{len=$3-$2; sum+=len} END {print "Average peak length after filtering: "sum/NR}' analysis/peaks2_mod/${sample}_broad_peaks_filtered.broadPeak

    # Blacklist filtering (optional)
    bedtools intersect -v \
        -a analysis/peaks2_mod/${sample}_broad_peaks_filtered.broadPeak \
        -b hg38-blacklist.v2.bed \
        > analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak

    # QC metrics (using final filtered peaks)
    if [ ! -s analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak ]; then
        echo "Error: Final filtered peaks file is empty!"
        exit 1
    fi

    total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
    reads_in_broad_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
        -b analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak -u | samtools view -c)
    broad_frip=$(echo "scale=4; $reads_in_broad_peaks / $total_reads" | bc)
    broad_peak_count=$(wc -l < analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak)
    broad_mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' \
        analysis/peaks2_mod/${sample}_broad_peaks_final.broadPeak)
    
    # Save QC metrics
    echo "${sample},broad,${total_reads},${reads_in_broad_peaks},${broad_frip},${broad_peak_count},${broad_mean_length}" \
        > analysis/qc2_mod/${sample}_broad_metrics.csv
fi
