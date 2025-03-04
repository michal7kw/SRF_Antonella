#!/usr/bin/env Rscript

# This script performs differential expression analysis using DESeq2
# It takes a count matrix and metadata file as input and produces DESeq2 results

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
})

# Get input and output file paths from snakemake
counts_file <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]
results_file <- snakemake@output[["results"]]
normalized_counts_file <- snakemake@output[["normalized_counts"]]
rds_file <- snakemake@output[["rds"]]
ma_plot_file <- snakemake@output[["ma_plot"]]
volcano_plot_file <- snakemake@output[["volcano_plot"]]
condition1 <- snakemake@params[["condition1"]]
condition2 <- snakemake@params[["condition2"]]
log_file <- snakemake@log[[1]]

# Set up logging
log <- file(log_file, open = "wt")
sink(log, type = "message")

message("Starting DESeq2 analysis for ", condition1, " vs ", condition2)
message("Reading count data from: ", counts_file)
message("Reading metadata from: ", metadata_file)

# Read the count matrix and metadata
count_matrix <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, check.names = FALSE)

# Filter metadata to include only the conditions we're comparing
metadata_filtered <- metadata %>%
  filter(condition %in% c(condition1, condition2))

# Filter count matrix to include only samples in the filtered metadata
count_matrix <- count_matrix[, metadata_filtered$sample]

# Ensure samples in count matrix and metadata are in the same order
metadata_filtered <- metadata_filtered %>%
  filter(sample %in% colnames(count_matrix)) %>%
  arrange(match(sample, colnames(count_matrix)))

# Verify that sample order matches
if (!all(metadata_filtered$sample == colnames(count_matrix))) {
  stop("Sample order in metadata does not match count matrix")
}

# Create DESeq2 dataset
message("Creating DESeq2 dataset")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata_filtered,
  design = ~ condition
)

# Set the reference level to condition2 (control)
dds$condition <- relevel(dds$condition, ref = condition2)

# Filter out genes with low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
message("Running DESeq2")
dds <- DESeq(dds)

# Get results
message("Extracting results")
res <- results(dds, contrast = c("condition", condition1, condition2))

# Add gene names to results
res$gene_id <- rownames(res)

# Convert to data frame and sort by adjusted p-value
res_df <- as.data.frame(res) %>%
  arrange(padj)

# Save normalized counts
message("Saving normalized counts")
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = normalized_counts_file, quote = FALSE)

# Save DESeq2 results
message("Saving DESeq2 results to: ", results_file)
write.csv(res_df, file = results_file, quote = FALSE)

# Save DESeq2 object for later use
message("Saving DESeq2 object to: ", rds_file)
saveRDS(dds, file = rds_file)

# Create MA plot
message("Creating MA plot")
pdf(ma_plot_file, width = 10, height = 8)
plotMA(res, main = paste0(condition1, " vs ", condition2), ylim = c(-5, 5))
dev.off()

# Create volcano plot
message("Creating volcano plot")
pdf(volcano_plot_file, width = 12, height = 10)
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0(condition1, " vs ", condition2),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                col = c('gray', 'blue', 'green', 'red'),
                colAlpha = 0.5,
                legendLabels = c('NS', 'Log2FC', 'p-value', 'p-value & Log2FC'),
                drawConnectors = TRUE,
                widthConnectors = 0.5)
dev.off()

message("DESeq2 analysis complete for ", condition1, " vs ", condition2)

# Close the log file
sink(type = "message")
close(log)
