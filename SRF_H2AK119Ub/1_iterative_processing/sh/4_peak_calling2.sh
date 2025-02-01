#!/bin/bash

#SBATCH --job-name=4_peak_calling2
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling/4_peak_calling2.err"
#SBATCH --output="logs/4_peak_calling/4_peak_calling2.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/{peaks2,qc2} 
#logs/peaks2

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
        --outdir analysis/peaks2 \
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

    check_output analysis/peaks2/${sample}_broad_peaks.broadPeak

    # Post-filtering (correct column 7 for fold-enrichment)
    awk 'BEGIN{OFS="\t"} {
        name=sprintf("%s_broad_peak_%d", "'${sample}'", NR);
        print $1, $2, $3, name, $5, $6, $7, $8, $9
    }' analysis/peaks2/${sample}_broad_peaks.broadPeak | \
    awk 'BEGIN{OFS="\t"} $7 >= 3 {print $0}' \
        > analysis/peaks2/${sample}_broad_peaks_filtered.broadPeak

    # Blacklist filtering (optional)
    bedtools intersect -v \
        -a analysis/peaks2/${sample}_broad_peaks_filtered.broadPeak \
        -b hg38-blacklist.v2.bed \
        > analysis/peaks2/${sample}_broad_peaks_final.broadPeak

    # QC metrics (using final filtered peaks)
    if [ ! -s analysis/peaks2/${sample}_broad_peaks_final.broadPeak ]; then
        echo "Error: Final filtered peaks file is empty!"
        exit 1
    fi

    total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
    reads_in_broad_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
        -b analysis/peaks2/${sample}_broad_peaks_final.broadPeak -u | samtools view -c)
    broad_frip=$(echo "scale=4; $reads_in_broad_peaks / $total_reads" | bc)
    broad_peak_count=$(wc -l < analysis/peaks2/${sample}_broad_peaks_final.broadPeak)
    broad_mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' \
        analysis/peaks2/${sample}_broad_peaks_final.broadPeak)
    
    # Save QC metrics
    echo "${sample},broad,${total_reads},${reads_in_broad_peaks},${broad_frip},${broad_peak_count},${broad_mean_length}" \
        > analysis/qc2/${sample}_broad_metrics.csv
fi

echo "Completed processing ${sample} for ${analysis_type} peak calling"
