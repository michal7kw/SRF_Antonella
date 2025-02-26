#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=6_bigwig          # Job name
#SBATCH --account=kubacki.michal      # Account to charge resources to
#SBATCH --mem=64GB                   # Memory (RAM) required per node
#SBATCH --time=INFINITE              # Maximum job runtime (unlimited)
#SBATCH --nodes=1                    # Number of nodes to use
#SBATCH --ntasks=32                  # Number of tasks (cores) to use
#SBATCH --array=0-5                  # Job array index (0 to 5, for 6 samples)
#SBATCH --mail-type=ALL              # Send email for all job events
#SBATCH --mail-user=kubacki.michal@hsr.it  # Email address for notifications
#SBATCH --error="logs/6_bigwig.err"  # Standard error log file
#SBATCH --output="logs/6_bigwig.out" # Standard output log file

# Documentation:
# This script generates a bigWig file and estimates library complexity for a given sample.
# It uses deeptools bamCoverage to create the bigWig file and Picard's EstimateLibraryComplexity
# to assess library complexity. The script is designed to be run as a SLURM job array.
# It processes each sample once, generating a bigWig file and a library complexity estimate.
#
# Input:
#   - analysis/3_alignment/${sample}.dedup.bam: Deduplicated BAM file for the sample.
#   - analysis/3_alignment/${sample}.dedup.bam.bai: Index file for the deduplicated BAM file.
#
# Output:
#   - analysis/6_bigwig/${sample}.bw: BigWig file containing normalized read coverage.
#   - analysis/6_bigwig/${sample}_complexity.txt: Text file containing library complexity metrics.
#   - logs/6_bigwig.err: Standard error log file for the SLURM job array.
#   - logs/6_bigwig.out: Standard output log file for the SLURM job array.

set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error.

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake  # Activate the 'snakemake' conda environment.

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR  # Change the current directory to the working directory.

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/6_bigwig"

# Create necessary directories if they don't exist.
mkdir -p ${OUTPUT_DIR} logs/peaks

# Define sample names (matching order with alignment script).
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)
# samples=(YAF_3) # This line is commented out, used for testing with a single sample.
# analysis_types=(broad narrow) # This line is commented out, used for different analysis types.
analysis_types=(broad) # Currently, only 'broad' analysis type is used.
sample_idx=$SLURM_ARRAY_TASK_ID  # Use the array task ID directly as the sample index
analysis_idx=0  # Always use 'broad' analysis type
sample=${samples[$sample_idx]}          # Get the sample name based on the SLURM array index.
analysis_type=${analysis_types[$analysis_idx]}  # Get the analysis type (currently unused, always 'broad').

# Function to check if input/output exists and is not empty.
check_output() {
    local file=$1  # The first argument is the file to check.
    if [[ ! -s $file ]]; then  # Check if the file does NOT exist or is empty (-s checks for non-zero size).
        echo "Error: Output file $file is empty or does not exist"
        exit 1  # Exit with an error code.
    fi
}

# Check if input BAM exists and is properly indexed.
if [[ ! -f analysis/3_alignment/${sample}.dedup.bam ]] || [[ ! -f analysis/3_alignment/${sample}.dedup.bam.bai ]]; then
    echo "Error: Input BAM or index not found for ${sample}"
    exit 1  # Exit with an error code if either the BAM file or its index is missing.
fi

echo "Processing ${sample} for ${analysis_type} peak calling..."

# Generate bigwig only once per sample (when analysis_idx is 0).
if [[ $analysis_idx -eq 0 ]]; then
    echo "Generating bigWig file for ${sample}..."
    # Use deeptools bamCoverage to generate a bigWig file.
    bamCoverage -b analysis/3_alignment/${sample}.dedup.bam \
        -o ${OUTPUT_DIR}/${sample}.bw \
        --binSize 10 \
        --normalizeUsing RPKM \
        --smoothLength 30 \
        --extendReads \
        --centerReads \
        --numberOfProcessors 32
    # -b: Input BAM file.
    # -o: Output bigWig file.
    # --binSize: Size of the bins (in base pairs) for averaging the coverage.
    # --normalizeUsing: Method for normalizing the coverage (RPKM: Reads Per Kilobase Million).
    # --smoothLength: Smoothing window size (in base pairs).
    # --extendReads: Extend reads to the fragment length.
    # --centerReads: Center reads before extending.
    # --numberOfProcessors: Number of processors to use.
    check_output ${OUTPUT_DIR}/${sample}.bw  # Check if the bigWig file was created successfully.

    echo "Calculating library complexity for ${sample}..."
    # Use Picard's EstimateLibraryComplexity to estimate library complexity.
    picard EstimateLibraryComplexity \
        I=analysis/3_alignment/${sample}.dedup.bam \
        O=${OUTPUT_DIR}/${sample}_complexity.txt
    # I: Input BAM file.
    # O: Output text file with complexity metrics.
    check_output ${OUTPUT_DIR}/${sample}_complexity.txt  # Check if the complexity file was created successfully.
fi

echo "Completed processing ${sample} for ${analysis_type}"