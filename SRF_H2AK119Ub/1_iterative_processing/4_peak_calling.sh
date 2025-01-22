#!/bin/bash
#SBATCH --job-name=4_peak_calling
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-11
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling_%a.err"
#SBATCH --output="logs/4_peak_calling_%a.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Changed to match alignment script

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/{peaks,visualization,qc} logs/peaks

# Define samples (matching order with alignment script)
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
analysis_types=(broad narrow)
sample_idx=$((SLURM_ARRAY_TASK_ID / 2))
analysis_idx=$((SLURM_ARRAY_TASK_ID % 2))
sample=${samples[$sample_idx]}
analysis_type=${analysis_types[$analysis_idx]}

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

# Generate bigwig only once per sample (when analysis_idx is 0)
if [[ $analysis_idx -eq 0 ]]; then
    echo "Generating bigWig file for ${sample}..."
    bamCoverage -b analysis/aligned/${sample}.dedup.bam \
        -o analysis/visualization/${sample}.bw \
        --binSize 10 \
        --normalizeUsing RPKM \
        --smoothLength 30 \
        --extendReads \
        --centerReads \
        --numberOfProcessors 32
    check_output analysis/visualization/${sample}.bw

    echo "Calculating library complexity for ${sample}..."
    picard EstimateLibraryComplexity \
        I=analysis/aligned/${sample}.dedup.bam \
        O=analysis/qc/${sample}_complexity.txt
    check_output analysis/qc/${sample}_complexity.txt
fi

# Peak calling based on analysis type
if [[ $analysis_type == "broad" ]]; then
    echo "Calling broad peaks for ${sample}..."
    macs2 callpeak \
        -t analysis/aligned/${sample}.dedup.bam \
        -f BAMPE \
        -g hs \
        --broad \
        -n ${sample}_broad \
        --outdir analysis/peaks \
        -q 0.05 \
        --broad-cutoff 0.1 \
        --keep-dup all \
        2> logs/peaks/${sample}_broad_peaks.log
    
    check_output analysis/peaks/${sample}_broad_peaks.broadPeak

    # Calculate broad peak QC metrics
    total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
    reads_in_broad_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
        -b analysis/peaks/${sample}_broad_peaks.broadPeak -u | samtools view -c)
    broad_frip=$(echo "scale=4; $reads_in_broad_peaks / $total_reads" | bc)
    broad_peak_count=$(wc -l < analysis/peaks/${sample}_broad_peaks.broadPeak)
    broad_mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' analysis/peaks/${sample}_broad_peaks.broadPeak)
    
    # Save broad peak QC metrics
    echo "${sample},broad,${total_reads},${reads_in_broad_peaks},${broad_frip},${broad_peak_count},${broad_mean_length}" > analysis/qc/${sample}_broad_metrics.csv
else
    echo "Calling narrow peaks for ${sample}..."
    macs2 callpeak \
        -t analysis/aligned/${sample}.dedup.bam \
        -f BAMPE \
        -g hs \
        -n ${sample}_narrow \
        --outdir analysis/peaks \
        -q 0.05 \
        --call-summits \
        --keep-dup all \
        2> logs/peaks/${sample}_narrow_peaks.log
    
    check_output analysis/peaks/${sample}_narrow_peaks.narrowPeak

    # Calculate narrow peak QC metrics
    total_reads=$(samtools view -c analysis/aligned/${sample}.dedup.bam)
    reads_in_narrow_peaks=$(bedtools intersect -a analysis/aligned/${sample}.dedup.bam \
        -b analysis/peaks/${sample}_narrow_peaks.narrowPeak -u | samtools view -c)
    narrow_frip=$(echo "scale=4; $reads_in_narrow_peaks / $total_reads" | bc)
    narrow_peak_count=$(wc -l < analysis/peaks/${sample}_narrow_peaks.narrowPeak)
    narrow_mean_length=$(awk '{sum += $3-$2} END {print sum/NR}' analysis/peaks/${sample}_narrow_peaks.narrowPeak)
    
    # Save narrow peak QC metrics
    echo "${sample},narrow,${total_reads},${reads_in_narrow_peaks},${narrow_frip},${narrow_peak_count},${narrow_mean_length}" > analysis/qc/${sample}_narrow_metrics.csv
fi

echo "Completed processing ${sample} for ${analysis_type} peak calling"