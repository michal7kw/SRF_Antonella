#!/bin/bash

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


PEAK_DIR="peaks2"

set -e  # Exit on error
set -u  # Exit on undefined variable

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR

# Define samples array
samples=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# Get current sample from array index
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Define input file
input_file="analysis/${PEAK_DIR}/${sample}_broad_peaks.broadPeak"
output_dir="analysis/peak_stats/${PEAK_DIR}"

mkdir -p ${output_dir}

# 1. Basic Statistics and Length Distribution
awk -v sample=$sample '
BEGIN {
    min_len=999999999;
    max_len=0;
    # Define length bins (adjust as needed)
    bin_size=500;
    max_bin=20000;
}
{
    len = $3 - $2;
    fold = $7;
    
    # Track min/max
    if(len < min_len) min_len = len;
    if(len > max_len) max_len = len;
    
    # Collect length distribution
    bin = int(len/bin_size) * bin_size;
    if(bin > max_bin) bin = max_bin;
    lengths[bin]++;
    
    # Collect fold enrichment stats
    sum_len += len;
    sum_fold += fold;
    
    # Store length and fold enrichment for correlation
    peaks[NR]["len"] = len;
    peaks[NR]["fold"] = fold;
    
    count++;
}
END {
    # Output summary statistics
    outfile = "'${output_dir}'/" sample "_peak_summary.txt"
    print "Sample: " sample > outfile
    print "Total Peaks: " count >> outfile
    print "Min Length: " min_len >> outfile
    print "Max Length: " max_len >> outfile
    print "Mean Length: " sum_len/count >> outfile
    print "Mean Fold Enrichment: " sum_fold/count >> outfile
    
    # Output length distribution
    lenfile = "'${output_dir}'/" sample "_length_dist.txt"
    for(bin in lengths) {
        print bin "\t" lengths[bin] > lenfile
    }
    
    # Output length vs fold enrichment
    foldfile = "'${output_dir}'/" sample "_length_vs_fold.txt"
    for(i=1; i<=count; i++) {
        print peaks[i]["len"] "\t" peaks[i]["fold"] > foldfile
    }
}' $input_file

# Sort the distribution file
sort -n -k1,1 ${output_dir}/${sample}_length_dist.txt > ${output_dir}/${sample}_length_dist_sorted.txt

# Create R script for visualization
cat << 'EOF' > ${output_dir}/plot_distributions.R
library(ggplot2)
library(gridExtra)

# Read data
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
data_dir <- args[2]

# Read length distribution
length_dist <- read.table(file.path(data_dir, paste0(sample, "_length_dist_sorted.txt")))
colnames(length_dist) <- c("Length", "Count")

# Read length vs fold enrichment
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

# Combine plots
pdf(file.path(data_dir, paste0(sample, "_peak_analysis.pdf")), width=12, height=8)
grid.arrange(p1, p2, ncol=2)
dev.off()

# Calculate and save quantiles
quantiles <- quantile(length_fold$Length, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
write.table(quantiles, file=file.path(data_dir, paste0(sample, "_length_quantiles.txt")),
            quote=FALSE, col.names=FALSE)
EOF

# Run R script
Rscript ${output_dir}/plot_distributions.R $sample ${output_dir}

# Display summary statistics
cat ${output_dir}/${sample}_peak_summary.txt
echo -e "\nLength Quantiles:"
cat ${output_dir}/${sample}_length_quantiles.txt 