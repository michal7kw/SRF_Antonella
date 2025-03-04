#!/usr/bin/env Rscript

# This script generates visualizations from DESeq2 results
# It creates PCA plots, sample correlation heatmaps, and heatmaps of top DEGs
# This is a standalone version that can be run independently of the Snakefile

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-d", "--dds_dir"), type="character", default="results/deseq2",
              help="Directory containing DESeq2 results [default: %default]"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
              help="Metadata file [default: %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="results/deseq2",
              help="Output directory for visualization files [default: %default]"),
  make_option(c("--pca_width"), type="numeric", default=10,
              help="Width of PCA plot in inches [default: %default]"),
  make_option(c("--pca_height"), type="numeric", default=8,
              help="Height of PCA plot in inches [default: %default]"),
  make_option(c("--heatmap_width"), type="numeric", default=10,
              help="Width of sample correlation heatmap in inches [default: %default]"),
  make_option(c("--heatmap_height"), type="numeric", default=8,
              help="Height of sample correlation heatmap in inches [default: %default]"),
  make_option(c("--top_degs_width"), type="numeric", default=12,
              help="Width of top DEGs heatmap in inches [default: %default]"),
  make_option(c("--top_degs_height"), type="numeric", default=14,
              help="Height of top DEGs heatmap in inches [default: %default]"),
  make_option(c("--top_n_degs"), type="numeric", default=50,
              help="Number of top DEGs to include in heatmap [default: %default]"),
  make_option(c("--pca_title"), type="character", default="PCA of RNA-seq Samples",
              help="Title for PCA plot [default: %default]"),
  make_option(c("--heatmap_title"), type="character", default="Sample Correlation Heatmap",
              help="Title for sample correlation heatmap [default: %default]"),
  make_option(c("--top_degs_title"), type="character", default="Top Differentially Expressed Genes",
              help="Title for top DEGs heatmap [default: %default]"),
  make_option(c("--pca_point_size"), type="numeric", default=3,
              help="Point size for PCA plot [default: %default]"),
  make_option(c("--pca_text_size"), type="numeric", default=4,
              help="Text size for PCA plot labels [default: %default]"),
  make_option(c("--pca_theme"), type="character", default="bw",
              help="Theme for PCA plot (bw, classic, minimal, etc.) [default: %default]"),
  make_option(c("--color_palette"), type="character", default="Set1",
              help="RColorBrewer palette for condition colors [default: %default]"),
  make_option(c("--verbose"), type="logical", default=TRUE,
              help="Print verbose output [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to print messages if verbose is TRUE
log_message <- function(message) {
  if (opt$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  }
}

log_message("Starting visualization generation")

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", opt$output_dir))
}

# Define output file paths
pca_plot_file <- file.path(opt$output_dir, "pca_plot.pdf")
heatmap_file <- file.path(opt$output_dir, "sample_correlation_heatmap.pdf")
top_degs_heatmap_file <- file.path(opt$output_dir, "top_degs_heatmap.pdf")

# Read metadata
log_message(paste("Reading metadata from:", opt$metadata))
metadata <- read.csv(opt$metadata, check.names = FALSE)

# Find all DESeq2 RDS files
log_message("Finding DESeq2 RDS files")
comparisons <- list.dirs(opt$dds_dir, full.names = TRUE, recursive = FALSE)
comparisons <- comparisons[grepl("_vs_", basename(comparisons))]
dds_files <- file.path(comparisons, "dds.rds")
dds_files <- dds_files[file.exists(dds_files)]

if (length(dds_files) == 0) {
  stop("No DESeq2 RDS files found in ", opt$dds_dir)
}

log_message(paste("Found", length(dds_files), "DESeq2 RDS files"))
for (file in dds_files) {
  log_message(paste("  -", file))
}

# Load the first DESeq2 object to get normalized counts for all samples
log_message(paste("Loading DESeq2 object from:", dds_files[1]))
dds <- readRDS(dds_files[1])

# Get variance stabilized transformed data for all samples
log_message("Performing variance stabilizing transformation")
vst <- vst(dds, blind = FALSE)

# Create PCA plot
log_message("Creating PCA plot")
pdf(pca_plot_file, width = opt$pca_width, height = opt$pca_height)
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = opt$pca_point_size) +
  geom_text_repel(size = opt$pca_text_size, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(opt$pca_title) +
  theme_get() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Apply selected theme
if (opt$pca_theme == "bw") {
  pca_plot <- pca_plot + theme_bw()
} else if (opt$pca_theme == "classic") {
  pca_plot <- pca_plot + theme_classic()
} else if (opt$pca_theme == "minimal") {
  pca_plot <- pca_plot + theme_minimal()
}

print(pca_plot)
dev.off()
log_message(paste("PCA plot saved to:", pca_plot_file))

# Create sample correlation heatmap
log_message("Creating sample correlation heatmap")
sample_dists <- dist(t(assay(vst)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst)
colnames(sample_dist_matrix) <- colnames(vst)

# Get condition colors
condition_colors <- brewer.pal(min(length(unique(metadata$condition)), 9), opt$color_palette)
names(condition_colors) <- unique(metadata$condition)
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- metadata$sample
ann_colors <- list(Condition = condition_colors)

pdf(heatmap_file, width = opt$heatmap_width, height = opt$heatmap_height)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = opt$heatmap_title,
         fontsize = 10)
dev.off()
log_message(paste("Sample correlation heatmap saved to:", heatmap_file))

# Collect top DEGs from all comparisons
log_message("Collecting top DEGs from all comparisons")
all_top_degs <- list()

for (dds_file in dds_files) {
  log_message(paste("Processing DESeq2 object:", dds_file))
  dds_obj <- readRDS(dds_file)
  
  # Get results
  res <- results(dds_obj)
  
  # Get top N DEGs by adjusted p-value
  top_degs <- rownames(res[order(res$padj), ])[1:opt$top_n_degs]
  all_top_degs <- c(all_top_degs, list(top_degs))
}

# Get unique top DEGs
unique_top_degs <- unique(unlist(all_top_degs))
log_message(paste("Number of unique top DEGs:", length(unique_top_degs)))

# Create heatmap of top DEGs
log_message("Creating heatmap of top DEGs")
# Extract normalized counts for top DEGs
top_degs_counts <- assay(vst)[unique_top_degs, ]

# Scale the counts for better visualization
top_degs_counts_scaled <- t(scale(t(top_degs_counts)))

pdf(top_degs_heatmap_file, width = opt$top_degs_width, height = opt$top_degs_height)
pheatmap(top_degs_counts_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         main = opt$top_degs_title,
         fontsize = 10)
dev.off()
log_message(paste("Top DEGs heatmap saved to:", top_degs_heatmap_file))

log_message("Visualization generation complete")
