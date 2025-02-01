#!/bin/bash

#SBATCH --job-name=4b_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/4_peak_calling/4b_bigwig%a.err"
#SBATCH --output="logs/4_peak_calling/4b_bigwig%a.out"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Changed to match alignment script

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Create necessary directories
mkdir -p analysis/visualization logs/peaks

# Define samples (matching order with alignment script)
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
# samples=(YAF_3)
# analysis_types=(broad narrow)
analysis_types=(broad)
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

echo "Completed processing ${sample} for ${analysis_type}"