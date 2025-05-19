#!/usr/bin/env Rscript

# Script to filter a differential_peaks.csv file for significant peaks based on FDR
# and optionally also by log fold change (LFC).
#
# This script:
# 1. Reads a CSV file containing differential peak analysis results.
# 2. Filters peaks based on a user-specified FDR threshold (output 1).
# 3. Optionally, further filters the FDR-significant peaks by an LFC threshold (output 2).
# 4. Writes the filtered peaks to new CSV files.

suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
})

# Define command line options
option_list <- list(
    make_option(c("-i", "--input_csv"), type="character", default=NULL,
                help="Path to the input differential_peaks.csv file (required)", metavar="character"),
    make_option(c("-o", "--output_csv_fdr"), type="character", default=NULL,
                help="Path for the output CSV file containing FDR significant peaks (required)", metavar="character"),
    make_option(c("-f", "--fdr_threshold"), type="double", default=0.05,
                help="FDR threshold for significance. Default is 0.05.", metavar="double"),
    make_option(c("--output_csv_lfc"), type="character", default=NULL,
                help="Path for the output CSV file containing FDR AND LFC significant peaks (optional). If provided, --lfc_threshold must also be set.", metavar="character"),
    make_option(c("-l", "--lfc_threshold"), type="double", default=1.0,
                help="Absolute log2 fold change threshold (e.g., 1 for a 2-fold change). Used if --output_csv_lfc is specified. Default is 1.0.", metavar="double"),
    make_option(c("--lfc_column_name"), type="character", default="Fold",
                help="Name of the log fold change column in the input CSV. Default is 'Fold'.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, description = "Filter differential peaks CSV for significant results based on FDR and optionally LFC.")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_csv) || is.null(opt$output_csv_fdr)){
    print_help(opt_parser)
    stop("Input CSV and at least the primary output CSV (output_csv_fdr) file paths must be supplied.", call.=FALSE)
}

# Validate LFC options if the second output file is requested
if (!is.null(opt$output_csv_lfc) && is.null(opt$lfc_threshold)) {
    print_help(opt_parser)
    stop("--lfc_threshold must be specified if --output_csv_lfc is used.", call.=FALSE)
}
if (is.null(opt$output_csv_lfc) && !is.null(opt$lfc_threshold) && opt$lfc_threshold != 1.0) {
    # Allow lfc_threshold to be set even if output_csv_lfc is not, but warn if it's not default
    # This scenario is not typical but doesn't break anything.
}


# Log message function
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

log_message("Starting significant peak filtering process.")
log_message(sprintf("Input CSV: %s", opt$input_csv))
log_message(sprintf("Output CSV for FDR significant peaks: %s", opt$output_csv_fdr))
log_message(sprintf("FDR threshold: %.4f", opt$fdr_threshold))
if (!is.null(opt$output_csv_lfc)) {
    log_message(sprintf("Output CSV for FDR & LFC significant peaks: %s", opt$output_csv_lfc))
    log_message(sprintf("Absolute LFC threshold: %.2f (using column '%s')", opt$lfc_threshold, opt$lfc_column_name))
}


# Read the differential peaks CSV file
log_message(sprintf("Reading differential peaks CSV from: %s", opt$input_csv))
tryCatch({
    all_peaks_df <- readr::read_csv(opt$input_csv, col_types = readr::cols(), show_col_types = FALSE)
}, error = function(e) {
    log_message(sprintf("Error reading CSV file: %s. Details: %s", opt$input_csv, e$message), level="ERROR")
    stop("Failed to read input CSV file.")
})

if (nrow(all_peaks_df) == 0) {
    log_message("Input CSV file is empty. No output files will be generated.", level="WARNING")
    file.create(opt$output_csv_fdr) # Create empty primary output
    log_message(sprintf("Empty output CSV (FDR) created at: %s", opt$output_csv_fdr))
    if (!is.null(opt$output_csv_lfc)) {
        file.create(opt$output_csv_lfc) # Create empty secondary output
        log_message(sprintf("Empty output CSV (FDR & LFC) created at: %s", opt$output_csv_lfc))
    }
    quit(save = "no", status = 0)
}

# Check if the 'FDR' column exists
if (!"FDR" %in% colnames(all_peaks_df)) {
    log_message("ERROR: 'FDR' column not found in the input CSV. Cannot filter for FDR significance.", level="ERROR")
    log_message(sprintf("Available columns: %s", paste(colnames(all_peaks_df), collapse=", ")), level="INFO")
    stop("Missing 'FDR' column in input CSV.")
}

# Filter for FDR significant peaks
log_message(sprintf("Filtering peaks with FDR < %.4f", opt$fdr_threshold))
fdr_significant_peaks_df <- all_peaks_df %>%
    dplyr::filter(FDR < opt$fdr_threshold)

num_total_peaks <- nrow(all_peaks_df)
num_fdr_significant_peaks <- nrow(fdr_significant_peaks_df)

log_message(sprintf("Found %d FDR significant peaks (FDR < %.4f) out of %d total peaks (%.2f%%).",
            num_fdr_significant_peaks, opt$fdr_threshold, num_total_peaks,
            ifelse(num_total_peaks > 0, (num_fdr_significant_peaks / num_total_peaks) * 100, 0)))

# Write the FDR significant peaks to the primary output CSV file
if (num_fdr_significant_peaks == 0) {
    log_message(sprintf("No significant peaks found at FDR < %.4f.", opt$fdr_threshold), level="INFO")
    if (num_total_peaks > 0) {
         readr::write_csv(fdr_significant_peaks_df, opt$output_csv_fdr) # Writes header even if 0 rows
    } else {
         file.create(opt$output_csv_fdr)
    }
    log_message(sprintf("Output CSV for FDR significant peaks (empty or header-only) created at: %s", opt$output_csv_fdr))
} else {
    log_message(sprintf("Writing %d FDR significant peaks to: %s", num_fdr_significant_peaks, opt$output_csv_fdr))
    tryCatch({
        readr::write_csv(fdr_significant_peaks_df, opt$output_csv_fdr)
    }, error = function(e) {
        log_message(sprintf("Error writing FDR significant peaks CSV: %s. Details: %s", opt$output_csv_fdr, e$message), level="ERROR")
        stop("Failed to write FDR significant peaks output CSV file.")
    })
}

# If a second output file for LFC filtering is requested
if (!is.null(opt$output_csv_lfc)) {
    if (num_fdr_significant_peaks == 0) {
        log_message("No FDR significant peaks to filter by LFC. Second output file will be empty or header-only.", level="INFO")
        if (num_total_peaks > 0 && "FDR" %in% colnames(all_peaks_df) && opt$lfc_column_name %in% colnames(all_peaks_df) ) {
            # Create an empty df with correct columns to write header
            empty_df_for_lfc <- fdr_significant_peaks_df %>% dplyr::filter(FALSE)
            readr::write_csv(empty_df_for_lfc, opt$output_csv_lfc)
        } else {
            file.create(opt$output_csv_lfc)
        }
        log_message(sprintf("Output CSV for FDR & LFC significant peaks (empty or header-only) created at: %s", opt$output_csv_lfc))
    } else {
        # Check if the LFC column exists
        if (!opt$lfc_column_name %in% colnames(fdr_significant_peaks_df)) {
            log_message(sprintf("ERROR: LFC column '%s' not found in the input CSV. Cannot filter by LFC.", opt$lfc_column_name), level="ERROR")
            log_message(sprintf("Available columns: %s", paste(colnames(fdr_significant_peaks_df), collapse=", ")), level="INFO")
            # Create an empty file for the LFC output as we cannot proceed
            file.create(opt$output_csv_lfc)
            log_message(sprintf("Empty output CSV (FDR & LFC) created at: %s due to missing LFC column.", opt$output_csv_lfc))
        } else {
            log_message(sprintf("Further filtering FDR significant peaks by absolute LFC (column '%s') >= %.2f", opt$lfc_column_name, opt$lfc_threshold))
            fdr_lfc_significant_peaks_df <- fdr_significant_peaks_df %>%
                dplyr::filter(abs(!!sym(opt$lfc_column_name)) >= opt$lfc_threshold)

            num_fdr_lfc_significant_peaks <- nrow(fdr_lfc_significant_peaks_df)

            log_message(sprintf("Found %d peaks significant by FDR AND LFC (abs(%s) >= %.2f) out of %d FDR significant peaks (%.2f%%).",
                        num_fdr_lfc_significant_peaks, opt$lfc_column_name, opt$lfc_threshold, num_fdr_significant_peaks,
                        ifelse(num_fdr_significant_peaks > 0, (num_fdr_lfc_significant_peaks / num_fdr_significant_peaks) * 100, 0)))

            if (num_fdr_lfc_significant_peaks == 0) {
                log_message(sprintf("No peaks found meeting both FDR < %.4f and abs(LFC) >= %.2f criteria.", opt$fdr_threshold, opt$lfc_threshold), level="INFO")
                 if (num_fdr_significant_peaks > 0) { # Check if there were FDR peaks to begin with
                    readr::write_csv(fdr_lfc_significant_peaks_df, opt$output_csv_lfc) # Writes header
                } else {
                    file.create(opt$output_csv_lfc)
                }
                log_message(sprintf("Output CSV for FDR & LFC significant peaks (empty or header-only) created at: %s", opt$output_csv_lfc))
            } else {
                log_message(sprintf("Writing %d FDR & LFC significant peaks to: %s", num_fdr_lfc_significant_peaks, opt$output_csv_lfc))
                tryCatch({
                    readr::write_csv(fdr_lfc_significant_peaks_df, opt$output_csv_lfc)
                }, error = function(e) {
                    log_message(sprintf("Error writing FDR & LFC significant peaks CSV: %s. Details: %s", opt$output_csv_lfc, e$message), level="ERROR")
                    stop("Failed to write FDR & LFC significant peaks output CSV file.")
                })
            }
        }
    }
}

log_message("Peak filtering script complete.")