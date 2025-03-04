#!/usr/bin/env Rscript

# This script prepares the count matrix from featureCounts output for DESeq2 analysis
# It reads the featureCounts output and creates a clean count matrix with gene IDs as rows and samples as columns

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# Get input and output file paths from snakemake
counts_file <- snakemake@input[["counts"]]
output_file <- snakemake@output[["matrix"]]
log_file <- snakemake@log[[1]]

# Set up logging
log <- file(log_file, open = "wt")
sink(log, type = "message")

message("Reading count data from: ", counts_file)

# Read the featureCounts output
counts_data <- read.table(counts_file, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)

# Extract only the gene ID and count columns
count_matrix <- counts_data %>%
  select(Geneid, contains("Aligned.sortedByCoord.out.bam")) %>%
  # Clean up column names to just have sample names
  rename_with(~ gsub(".*star/(.+)/.+_Aligned.sortedByCoord.out.bam", "\\1", .), -Geneid)

# Set gene IDs as row names
rownames(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL

# Write the clean count matrix to a CSV file
message("Writing count matrix to: ", output_file)
write.csv(count_matrix, file = output_file, quote = FALSE)

message("Count matrix preparation complete")

# Close the log file
sink(type = "message")
close(log)
