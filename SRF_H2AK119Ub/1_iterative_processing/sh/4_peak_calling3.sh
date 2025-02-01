#!/bin/bash

#SBATCH --job-name=4_peak_calling3
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling/4_peak_calling3.err"
#SBATCH --output="logs/4_peak_calling/4_peak_calling3.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/{peaks3,qc3} 
#logs/peaks3

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
    echo "Calling broad peaks3 for ${sample} with stringent parameters..."
    macs2 callpeak \
        -t analysis/aligned/${sample}.dedup.bam \
        -f BAMPE \
        -g hs \
        --broad \
        -n ${sample}_broad \
        --outdir analysis/peaks3 \
        -q 0.05 \
        --broad-cutoff 0.05 \
        --keep-dup 1 \
        --min-length 500 \
        --max-gap 150 \
        --bdg \
        --buffer-size 100000

    check_output analysis/peaks3/${sample}_broad_peaks.broadPeak

    # Post-filtering of peaks
    # Filter by fold-enrichment (column 5)
    awk 'BEGIN{OFS="\t"} $5 >= 10 {print $0}' \
        analysis/peaks3/${sample}_broad_peaks.broadPeak > \
        analysis/peaks3/${sample}_broad_peaks_filtered.broadPeak

    # Add error checking for filtered peaks file
    check_output analysis/peaks3/${sample}_broad_peaks_filtered.broadPeak

    # Calculate QC metrics using filtered peaks
    total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
    reads_in_broad_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
        -b analysis/peaks3/${sample}_broad_peaks_filtered.broadPeak -u | samtools view -c)
    broad_frip=$(echo "scale=4; $reads_in_broad_peaks / $total_reads" | bc)
    broad_peak_count=$(wc -l < analysis/peaks3/${sample}_broad_peaks_filtered.broadPeak)
    broad_mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' \
        analysis/peaks3/${sample}_broad_peaks_filtered.broadPeak)
    
    # Save QC metrics
    echo "${sample},broad,${total_reads},${reads_in_broad_peaks},${broad_frip},${broad_peak_count},${broad_mean_length}" \
        > analysis/qc3/${sample}_broad_metrics.csv
fi

