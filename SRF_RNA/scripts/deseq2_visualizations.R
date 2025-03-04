#!/usr/bin/env Rscript

# This script generates additional visualizations from DESeq2 results
# It creates PCA plots, sample correlation heatmaps, and heatmaps of top DEGs

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Get input and output file paths from snakemake
dds_files <- snakemake@input[["dds_files"]]
metadata_file <- snakemake@input[["metadata"]]
pca_plot_file <- snakemake@output[["pca"]]
heatmap_file <- snakemake@output[["heatmap"]]
top_degs_heatmap_file <- snakemake@output[["top_degs_heatmap"]]
log_file <- snakemake@log[[1]]

# Set up logging
log <- file(log_file, open = "wt")
sink(log, type = "message")

message("Starting visualization generation")
message("Reading metadata from: ", metadata_file)

# Read metadata
metadata <- read.csv(metadata_file, check.names = FALSE)

# Load the first DESeq2 object to get normalized counts for all samples
message("Loading DESeq2 object from: ", dds_files[1])
dds <- readRDS(dds_files[1])

# Get variance stabilized transformed data for all samples
message("Performing variance stabilizing transformation")
vst <- vst(dds, blind = FALSE)

# Create PCA plot
message("Creating PCA plot")
pdf(pca_plot_file, width = 10, height = 8)
plotPCA(vst, intgroup = "condition") +
  theme_bw() +
  ggtitle("PCA of RNA-seq Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
dev.off()

# Create sample correlation heatmap
message("Creating sample correlation heatmap")
sample_dists <- dist(t(assay(vst)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst)
colnames(sample_dist_matrix) <- colnames(vst)

# Get condition colors
condition_colors <- brewer.pal(length(unique(metadata$condition)), "Set1")
names(condition_colors) <- unique(metadata$condition)
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- metadata$sample
ann_colors <- list(Condition = condition_colors)

pdf(heatmap_file, width = 10, height = 8)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Sample Correlation Heatmap",
         fontsize = 10)
dev.off()

# Collect top DEGs from all comparisons
message("Collecting top DEGs from all comparisons")
all_top_degs <- list()

for (dds_file in dds_files) {
  message("Processing DESeq2 object: ", dds_file)
  dds_obj <- readRDS(dds_file)
  
  # Get results
  res <- results(dds_obj)
  
  # Get top 50 DEGs by adjusted p-value
  top_degs <- rownames(res[order(res$padj), ])[1:50]
  all_top_degs <- c(all_top_degs, list(top_degs))
}

# Get unique top DEGs
unique_top_degs <- unique(unlist(all_top_degs))
message("Number of unique top DEGs: ", length(unique_top_degs))

# Create heatmap of top DEGs
message("Creating heatmap of top DEGs")
# Extract normalized counts for top DEGs
top_degs_counts <- assay(vst)[unique_top_degs, ]

# Scale the counts for better visualization
top_degs_counts_scaled <- t(scale(t(top_degs_counts)))

pdf(top_degs_heatmap_file, width = 12, height = 14)
pheatmap(top_degs_counts_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         main = "Top Differentially Expressed Genes",
         fontsize = 10)
dev.off()

message("Visualization generation complete")

# Close the log file
sink(type = "message")
close(log)
