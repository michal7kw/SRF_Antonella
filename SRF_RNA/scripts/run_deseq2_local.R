#!/usr/bin/env Rscript

# This script performs differential expression analysis using DESeq2
# It takes a count matrix and metadata file as input and produces DESeq2 results
# This is a standalone version adapted for local execution.

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
})

# Set working directory (assuming the script is in SRF_RNA/scripts)
# Absolute path to the project root: /mnt/d/Github/SRF_H2AK119Ub_cross_V5
# Script location: /mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_RNA/scripts
project_root <- "/mnt/d/Github/SRF_H2AK119Ub_cross_V5"
script_dir_relative_to_project_root <- "SRF_RNA/scripts"
full_script_dir <- file.path(project_root, script_dir_relative_to_project_root)

# Check if the script is being sourced or run directly, and set wd accordingly
if (interactive() || sys.nframe() == 0) { # sys.nframe() == 0 means script is run with Rscript
    # If run directly or sourced in interactive session, set working directory
    if (getwd() != full_script_dir) {
        setwd(full_script_dir)
        cat("Working directory set to:", getwd(), "\n")
    } else {
        cat("Working directory is already:", getwd(), "\n")
    }
}


opt <- list(
  counts_file = "../results/counts/all_samples_counts.txt", # Relative to SRF_RNA/scripts
  metadata_file = "../metadata.csv",                       # Relative to SRF_RNA/scripts
  output_base_dir = "../results/deseq2_local",             # Base directory for outputs, relative
  comparison_name = "YAF_vs_GFP",                          # Subdirectory for this specific comparison
  condition1 = "YAF",
  condition2 = "GFP",
  verbose = TRUE,
  ma_plot_width = 10,
  ma_plot_height = 8,
  volcano_plot_width = 12,
  volcano_plot_height = 10
)

# Construct full output directory path
opt$output_dir <- file.path(opt$output_base_dir, opt$comparison_name)

# Define output file paths based on opt$output_dir
opt$results_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_deseq2_results.csv"))
opt$normalized_counts_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_normalized_counts.csv"))
opt$rds_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_dds.rds"))
opt$ma_plot_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_ma_plot.pdf"))
opt$volcano_plot_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_volcano_plot.pdf"))
# opt$log_file <- file.path(opt$output_dir, paste0(opt$comparison_name, "_run_deseq2.log")) # Log file path if needed

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  cat(paste("Created output directory:", opt$output_dir, "\n"))
}

# Function to print messages if verbose is TRUE
log_message <- function(message_text) {
  if (opt$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message_text, "\n"))
  }
}

log_message(paste0("Starting DESeq2 analysis for ", opt$condition1, " vs ", opt$condition2))
log_message(paste0("Reading count data from: ", opt$counts_file))
log_message(paste0("Reading metadata from: ", opt$metadata_file))

# Read the count matrix and metadata
# Ensure paths are correct relative to the new working directory (SRF_RNA/scripts)
count_matrix_path <- opt$counts_file # Already relative
metadata_path <- opt$metadata_file   # Already relative

count_matrix <- read.delim(count_matrix_path, row.names = 1, check.names = FALSE, comment.char = "#")

# Clean column names of count_matrix
# Original colnames are like: results/star/C1/C1_Aligned.sortedByCoord.out.bam
# We need to extract the sample ID part, e.g., "C1"
if (ncol(count_matrix) > 0) { # Ensure there are columns to process
  original_colnames <- colnames(count_matrix)
  # Extract the directory name that is the sample ID
  # e.g., basename(dirname("results/star/C1/C1_Aligned.sortedByCoord.out.bam")) -> "C1"
  cleaned_colnames <- basename(dirname(original_colnames))
  colnames(count_matrix) <- cleaned_colnames
  log_message(paste0("Cleaned count matrix column names. Example: ", paste(head(cleaned_colnames), collapse=", ")))
} else {
  log_message("Warning: Count matrix has no columns after reading. Check file format and path.")
}
metadata <- read.csv(metadata_path, check.names = FALSE)

# Filter metadata to include only the conditions we're comparing
metadata_filtered <- metadata %>%
  filter(condition %in% c(opt$condition1, opt$condition2))

# Filter count matrix to include only samples in the filtered metadata
# Ensure sample column in metadata is named 'sample' or adjust as needed.
# Assuming metadata has a 'sample' column that matches colnames(count_matrix)
if (!"sample" %in% colnames(metadata_filtered)) {
    stop("Metadata must contain a 'sample' column matching count matrix column names.")
}
count_matrix_filtered <- count_matrix[, metadata_filtered$sample]


# Ensure samples in count matrix and metadata are in the same order
metadata_filtered <- metadata_filtered %>%
  filter(sample %in% colnames(count_matrix_filtered)) %>%
  arrange(match(sample, colnames(count_matrix_filtered)))

# Verify that sample order matches
if (!all(metadata_filtered$sample == colnames(count_matrix_filtered))) {
  stop("Sample order in metadata does not match count matrix after filtering.")
}

# Create DESeq2 dataset
log_message("Creating DESeq2 dataset")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered, # Use filtered count matrix
  colData = metadata_filtered,
  design = ~ condition
)

# Set the reference level to condition2 (control)
dds$condition <- relevel(dds$condition, ref = opt$condition2)

# Filter out genes with low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
log_message("Running DESeq2")
dds <- DESeq(dds)

# Get results
log_message("Extracting results")
res <- results(dds, contrast = c("condition", opt$condition1, opt$condition2))

# Add gene names to results
res$gene_id <- rownames(res)

# Convert to data frame and sort by adjusted p-value
res_df <- as.data.frame(res) %>%
  arrange(padj)

# Save normalized counts
log_message(paste0("Saving normalized counts to: ", opt$normalized_counts_file))
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = opt$normalized_counts_file, quote = FALSE)

# Save DESeq2 results
log_message(paste0("Saving DESeq2 results to: ", opt$results_file))
write.csv(res_df, file = opt$results_file, quote = FALSE)

# Save DESeq2 object for later use
log_message(paste0("Saving DESeq2 object to: ", opt$rds_file))
saveRDS(dds, file = opt$rds_file)

# Create MA plot
log_message(paste0("Creating MA plot: ", opt$ma_plot_file))
pdf(opt$ma_plot_file, width = opt$ma_plot_width, height = opt$ma_plot_height)
plotMA(res, main = paste0(opt$condition1, " vs ", opt$condition2), ylim = c(-5, 5))
dev.off()

# Create volcano plot
log_message(paste0("Creating volcano plot: ", opt$volcano_plot_file))
pdf(opt$volcano_plot_file, width = opt$volcano_plot_width, height = opt$volcano_plot_height)
EnhancedVolcano(res_df,
                lab = res_df$gene_id, # Corrected from rownames(res_df)
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0(opt$condition1, " vs ", opt$condition2),
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

log_message(paste0("DESeq2 analysis complete for ", opt$condition1, " vs ", opt$condition2))
log_message(paste0("Output files are in: ", opt$output_dir))