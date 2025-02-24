#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=analyze_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/analyze_peaks.err"
#SBATCH --output="logs/analyze_peaks.out"
#SBATCH --array=0-5

# Documentation:
# This script performs comprehensive analysis of ChIP-seq peak files. It includes:
#   - Calculation of basic peak statistics (total count, minimum/maximum/mean length, mean fold enrichment).
#   - Analysis of peak length distribution.
#   - Correlation analysis between peak length and fold enrichment.
#   - Visualization of peak length distribution and length vs. fold enrichment.
#
# Inputs:
#   - BroadPeak files: Located in analysis/5_peak_calling/, named as {sample}_broad_peaks_final.broadPeak.
#     These files should contain ChIP-seq peak data in BroadPeak format.
#   - Sample names: Defined in the 'samples' array (GFP_1, GFP_2, GFP_3, YAF_1, YAF_2, YAF_3). The script uses the SLURM_ARRAY_TASK_ID
#     to process each sample in parallel.
#
# Outputs:
#   - Text files (in analysis/13_analyze_peaks/):
#     - peak_summary.txt: Contains summary statistics for each sample (total peaks, min/max/mean length, mean fold enrichment).
#     - length_dist.txt:  Contains the distribution of peak lengths in 500bp bins.
#     - length_vs_fold.txt: Contains peak length and fold enrichment data for each peak, used for correlation analysis.
#     - length_dist_sorted.txt: Sorted version of length_dist.txt, used for plotting.
#     - length_quantiles.txt: Contains the 5th, 25th, 50th, 75th, and 95th percentiles of peak lengths.
#   - PDF files (in analysis/13_analyze_peaks/):
#     - peak_analysis.pdf: A PDF containing two plots: peak length distribution (histogram) and length vs fold enrichment (scatter plot).
#   - Directory: analysis/13_analyze_peaks/. This is the main output directory where all result files are stored.

# Directory containing peak files
PEAK_DIR="5_peak_calling"

# Error handling
set -e  # Exit on error
set -u  # Exit on undefined variable

# Activate conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define and navigate to working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Define sample names array
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# Get current sample from SLURM array index
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Define the main output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/13_analyze_peaks"

# Define input and output paths
input_file="analysis/${PEAK_DIR}/${sample}_broad_peaks_final.broadPeak"
output_dir="${OUTPUT_DIR}"

# Create output directory if it doesn't exist
mkdir -p ${output_dir}

# Step 1: Basic Statistics and Length Distribution
# Process BroadPeak file to calculate:
# - Total number of peaks
# - Minimum, maximum, and mean peak length
# - Mean fold enrichment
# - Length distribution in 500bp bins
# - Length vs fold enrichment correlation
awk -v sample=$sample -v output_dir="${output_dir}" '
BEGIN {
    # Initialize variables
    min_len=999999999;
    max_len=0;
    bin_size=500;       # Bin size for length distribution
    max_bin=20000;      # Maximum bin value
}
{
    # Calculate peak length and get fold enrichment
    len = $3 - $2;
    fold = $7;
    
    # Track minimum and maximum lengths
    if(len < min_len) min_len = len;
    if(len > max_len) max_len = len;
    
    # Bin lengths for distribution analysis
    bin = int(len/bin_size) * bin_size;
    if(bin > max_bin) bin = max_bin;
    lengths[bin]++;
    
    # Accumulate values for mean calculations
    sum_len += len;
    sum_fold += fold;
    
    # Store individual peak data for correlation analysis
    peaks[NR]["len"] = len;
    peaks[NR]["fold"] = fold;
    
    count++;
}
END {
    # Output summary statistics
    outfile = output_dir "/" sample "_peak_summary.txt"
    print "Sample: " sample > outfile
    print "Total Peaks: " count >> outfile
    print "Min Length: " min_len >> outfile
    print "Max Length: " max_len >> outfile
    print "Mean Length: " sum_len/count >> outfile
    print "Mean Fold Enrichment: " sum_fold/count >> outfile
    
    # Output length distribution data
    lenfile = output_dir "/" sample "_length_dist.txt"
    for(bin in lengths) {
        print bin "\t" lengths[bin] > lenfile
    }
    
    # Output length vs fold enrichment data
    foldfile = output_dir "/" sample "_length_vs_fold.txt"
    for(i=1; i<=count; i++) {
        print peaks[i]["len"] "\t" peaks[i]["fold"] > foldfile
    }
}' $input_file

# Sort length distribution file for plotting
sort -n -k1,1 ${output_dir}/${sample}_length_dist.txt > ${output_dir}/${sample}_length_dist_sorted.txt

# Step 2: Create R script for visualization
# Generates two plots:
# 1. Peak length distribution (histogram)
# 2. Length vs fold enrichment (scatter plot)
cat << EOF > ${output_dir}/plot_distributions.R
library(ggplot2)
library(gridExtra)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
data_dir <- args[2]

# Read and process length distribution data
length_dist <- read.table(file.path(data_dir, paste0(sample, "_length_dist_sorted.txt")))
colnames(length_dist) <- c("Length", "Count")

# Read and process length vs fold enrichment data
length_fold <- read.table(file.path(data_dir, paste0(sample, "_length_vs_fold.txt")))
colnames(length_fold) <- c("Length", "FoldEnrichment")

# Create length distribution plot
p1 <- ggplot(length_dist, aes(x=Length, y=Count)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(title=paste("Peak Length Distribution -", sample),
         x="Peak Length (bp)", y="Number of Peaks") +
    scale_x_continuous(breaks=seq(0, max(length_dist$Length), by=2000))

# Create length vs fold enrichment scatter plot
p2 <- ggplot(length_fold, aes(x=Length, y=FoldEnrichment)) +
    geom_point(alpha=0.3) +
    theme_minimal() +
    labs(title=paste("Peak Length vs Fold Enrichment -", sample),
         x="Peak Length (bp)", y="Fold Enrichment") +
    scale_x_continuous(breaks=seq(0, max(length_fold$Length), by=2000))

# Combine plots into single PDF
pdf(file.path(data_dir, paste0(sample, "_peak_analysis.pdf")), width=12, height=8)
grid.arrange(p1, p2, ncol=2)
dev.off()

# Calculate and save length quantiles
quantiles <- quantile(length_fold$Length, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
write.table(quantiles, file=file.path(data_dir, paste0(sample, "_length_quantiles.txt")),
            quote=FALSE, col.names=FALSE)
EOF

# Step 3: Run R script to generate visualizations
Rscript ${output_dir}/plot_distributions.R $sample ${output_dir}

# Step 4: Display summary statistics
cat ${output_dir}/${sample}_peak_summary.txt
echo -e "\nLength Quantiles:"
cat ${output_dir}/${sample}_length_quantiles.txt