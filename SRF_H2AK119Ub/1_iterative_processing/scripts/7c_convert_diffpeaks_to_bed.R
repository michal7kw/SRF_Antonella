#!/usr/bin/env Rscript

# Script to convert differential_peaks.csv from DiffBind to BED format for IGV visualization
#
# This script reads a CSV file containing differential peak analysis results (typically from DiffBind),
# extracts relevant columns, formats them into BED6 (chrom, chromStart, chromEnd, name, score, strand),
# and writes the output to a BED file.
#
# BED format details:
# - chrom: Chromosome name
# - chromStart: 0-based start position
# - chromEnd: 1-based end position
# - name: A name for the feature
# - score: A score (e.g., fold change, p-value)
# - strand: +, -, or . (for not stranded)

# Load necessary libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
})

# Define command line options
option_list <- list(
    make_option(c("-i", "--input_csv"), type="character", default=NULL,
                help="Path to the differential_peaks.csv file (required)", metavar="character"),
    make_option(c("-o", "--output_bed"), type="character", default=NULL,
                help="Path for the output BED file (required)", metavar="character"),
    make_option(c("-s", "--score_column"), type="character", default="Fold",
                help="Column name from CSV to use for BED score (e.g., Fold, FDR). Default is 'Fold'.", metavar="character"),
    make_option(c("--add_chr_prefix"), type="logical", default=FALSE, action="store_true",
                help="Add 'chr' prefix to chromosome names if not present (for UCSC style genome browsers). Default is FALSE.")
)

opt_parser <- OptionParser(option_list=option_list, description = "Convert DiffBind CSV output to BED format for IGV.")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_csv) || is.null(opt$output_bed)){
    print_help(opt_parser)
    stop("Input CSV and output BED file paths must be supplied.", call.=FALSE)
}

# Log message function
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

log_message(sprintf("Starting BED conversion process."))
log_message(sprintf("Input CSV: %s", opt$input_csv))
log_message(sprintf("Output BED: %s", opt$output_bed))
log_message(sprintf("Score column: %s", opt$score_column))
log_message(sprintf("Add 'chr' prefix: %s", opt$add_chr_prefix))

# Read the differential peaks CSV file
log_message(sprintf("Reading differential peaks CSV from: %s", opt$input_csv))
tryCatch({
    diff_peaks_df <- readr::read_csv(opt$input_csv, col_types = readr::cols())
}, error = function(e) {
    log_message(sprintf("Error reading CSV file: %s. Details: %s", opt$input_csv, e$message), level="ERROR")
    stop("Failed to read input CSV file.")
})

if (nrow(diff_peaks_df) == 0) {
    log_message("Input CSV file is empty. No BED file will be generated.", level="WARNING")
    quit(save = "no", status = 0)
}

# Validate that the specified score column exists
if (!opt$score_column %in% colnames(diff_peaks_df)) {
    available_cols <- paste(colnames(diff_peaks_df), collapse=", ")
    log_message(sprintf("Specified score column '%s' not found in CSV. Available columns: [%s]",
                opt$score_column, available_cols), level="ERROR")
    stop("Invalid score column specified. Please check the column names in your CSV file.")
}

# Ensure essential columns from DiffBind output are present
required_cols <- c("seqnames", "start", "end", "strand")
missing_cols <- setdiff(required_cols, colnames(diff_peaks_df))
if (length(missing_cols) > 0) {
    log_message(sprintf("Essential columns missing from CSV: %s. These are needed for BED conversion.",
                paste(missing_cols, collapse=", ")), level="ERROR")
    stop("Missing essential columns in input CSV.")
}


log_message("Formatting data into BED6 format (chrom, chromStart, chromEnd, name, score, strand).")

# BED format: chrom, chromStart (0-based), chromEnd (1-based), name, score, strand
# DiffBind output (GRanges to data.frame): seqnames, start (1-based), end (1-based), strand, Fold, FDR, etc.

bed_df <- diff_peaks_df %>%
    dplyr::mutate(
        chrom = as.character(seqnames),
        chromStart = as.integer(start) - 1, # Convert to 0-based start
        chromEnd = as.integer(end),
        # Create a unique name for each feature, incorporating original coordinates and score
        name = paste0(chrom, ":", start, "-", end, "_", opt$score_column, ":", round(as.numeric(.data[[opt$score_column]]), 2)),
        score = as.numeric(.data[[opt$score_column]]), # Use user-specified score column
        strand = ifelse(strand == "*" | is.na(strand), ".", as.character(strand)) # Ensure strand is ., +, or -
    ) %>%
    dplyr::filter(!is.na(chrom) & !is.na(chromStart) & !is.na(chromEnd)) # Remove rows with NA in critical BED columns

if (nrow(bed_df) == 0) {
    log_message("No valid rows remaining after initial filtering (e.g., NA in chr/start/end). BED file will be empty or not generated.", level="WARNING")
    # Optionally, write an empty BED file or just exit
    file.create(opt$output_bed) # Create an empty file
    log_message(sprintf("Empty BED file created at: %s", opt$output_bed))
    quit(save = "no", status = 0)
}


# Add "chr" prefix if requested and if not already present
if (opt$add_chr_prefix) {
    log_message("Adding 'chr' prefix to chromosome names where needed.")
    bed_df <- bed_df %>%
        dplyr::mutate(chrom = ifelse(!startsWith(chrom, "chr"), paste0("chr", chrom), chrom))
}

# Select and order columns for BED6 format
bed_df_final <- bed_df %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand)

# Ensure chromStart is not negative (can happen if original start was 0 or 1 and not handled)
if(any(bed_df_final$chromStart < 0)) {
    n_negative <- sum(bed_df_final$chromStart < 0)
    log_message(sprintf("%d records found with negative chromStart values. Clamping to 0.", n_negative), level="WARNING")
    bed_df_final <- bed_df_final %>% dplyr::mutate(chromStart = ifelse(chromStart < 0, 0, chromStart))
}

# Ensure chromStart < chromEnd
invalid_intervals <- bed_df_final %>% dplyr::filter(chromStart >= chromEnd)
if (nrow(invalid_intervals) > 0) {
    log_message(sprintf("%d records found where chromStart >= chromEnd. These will be removed.", nrow(invalid_intervals)), level="WARNING")
    bed_df_final <- bed_df_final %>% dplyr::filter(chromStart < chromEnd)
}

if (nrow(bed_df_final) == 0) {
    log_message("No valid rows remaining after all filtering steps. BED file will be empty or not generated.", level="WARNING")
    file.create(opt$output_bed) # Create an empty file
    log_message(sprintf("Empty BED file created at: %s", opt$output_bed))
    quit(save = "no", status = 0)
}

# Sort the BED data by chromosome then start position
log_message("Sorting BED file by chromosome and start position.")
bed_df_final <- bed_df_final %>%
    dplyr::arrange(chrom, chromStart)

log_message(sprintf("Writing sorted BED file with %d records to: %s", nrow(bed_df_final), opt$output_bed))
tryCatch({
    readr::write_tsv(bed_df_final,
                     file = opt$output_bed,
                     col_names = FALSE)
}, error = function(e) {
    log_message(sprintf("Error writing BED file: %s. Details: %s", opt$output_bed, e$message), level="ERROR")
    stop("Failed to write output BED file.")
})

log_message("BED file generation complete.")
log_message(sprintf("To load into IGV: File > Load from File... > Select '%s'", opt$output_bed))
