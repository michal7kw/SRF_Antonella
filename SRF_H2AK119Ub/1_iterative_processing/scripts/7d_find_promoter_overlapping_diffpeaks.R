#!/usr/bin/env Rscript

# Script to identify differential peaks overlapping with gene promoter regions
#
# This script:
# 1. Reads differential peak data (CSV from DiffBind).
# 2. Reads a gene annotation GTF file.
# 3. Defines promoter regions based on Transcription Start Sites (TSS) from the GTF.
# 4. Finds differential peaks that overlap these promoter regions.
# 5. Outputs the overlapping peaks in BED format.

suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomeInfoDb) # For seqlevelsStyle
})

# Define command line options
option_list <- list(
    make_option(c("-d", "--diff_peaks_csv"), type="character", default=NULL,
                help="Path to the differential_peaks.csv file (required)", metavar="character"),
    make_option(c("-g", "--gtf_file"), type="character", default=NULL,
                help="Path to the gene annotation GTF file (required)", metavar="character"),
    make_option(c("-o", "--output_bed"), type="character", default=NULL,
                help="Path for the output BED file of promoter-overlapping peaks (required)", metavar="character"),
    make_option(c("--promoter_upstream"), type="integer", default=2000,
                help="Distance upstream of TSS to define promoter region (bp). Default is 2000.", metavar="integer"),
    make_option(c("--promoter_downstream"), type="integer", default=500,
                help="Distance downstream of TSS to define promoter region (bp). Default is 500.", metavar="integer"),
    make_option(c("-s", "--score_column"), type="character", default="Fold",
                help="Column name from CSV to use for BED score (e.g., Fold, FDR). Default is 'Fold'.", metavar="character"),
    make_option(c("--add_chr_prefix_output"), type="logical", default=FALSE, action="store_true",
                help="Add 'chr' prefix to chromosome names in the output BED file. Default is FALSE."),
    make_option(c("--peak_chr_style"), type="character", default="NCBI",
                help="Chromosome style of input peaks ('NCBI' for 1,2,X or 'UCSC' for chr1,chr2,chrX). Default NCBI.", metavar="character"),
    make_option(c("--gtf_chr_style"), type="character", default="UCSC",
                help="Chromosome style of GTF file ('NCBI' or 'UCSC'). Default UCSC.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, description = "Find differential peaks overlapping gene promoters.")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$diff_peaks_csv) || is.null(opt$gtf_file) || is.null(opt$output_bed)){
    print_help(opt_parser)
    stop("Input CSV, GTF file, and output BED file paths must be supplied.", call.=FALSE)
}

# Log message function
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

log_message("Starting promoter overlap analysis.")
log_message(sprintf("Differential peaks CSV: %s", opt$diff_peaks_csv))
log_message(sprintf("GTF file: %s", opt$gtf_file))
log_message(sprintf("Output BED: %s", opt$output_bed))
log_message(sprintf("Promoter region: %d bp upstream, %d bp downstream of TSS", opt$promoter_upstream, opt$promoter_downstream))
log_message(sprintf("Score column for BED: %s", opt$score_column))
log_message(sprintf("Add 'chr' prefix to output BED: %s", opt$add_chr_prefix_output))
log_message(sprintf("Input peak chr style: %s", opt$peak_chr_style))
log_message(sprintf("Input GTF chr style: %s", opt$gtf_chr_style))


# 1. Read differential peaks CSV
log_message(sprintf("Reading differential peaks CSV from: %s", opt$diff_peaks_csv))
tryCatch({
    diff_peaks_df <- readr::read_csv(opt$diff_peaks_csv, col_types = readr::cols(), show_col_types = FALSE)
    # Ensure required columns exist
    required_peak_cols <- c("seqnames", "start", "end", "strand", opt$score_column)
    missing_peak_cols <- setdiff(required_peak_cols, colnames(diff_peaks_df))
    if (length(missing_peak_cols) > 0) {
        stop(paste("Missing required columns in peaks CSV:", paste(missing_peak_cols, collapse=", ")))
    }
    # Convert to GRanges
    peaks_gr <- GRanges(
        seqnames = Rle(diff_peaks_df$seqnames),
        ranges = IRanges(start = as.integer(diff_peaks_df$start), end = as.integer(diff_peaks_df$end)),
        strand = Rle(diff_peaks_df$strand),
        score = as.numeric(diff_peaks_df[[opt$score_column]]),
        peak_id = paste0(diff_peaks_df$seqnames, ":", diff_peaks_df$start, "-", diff_peaks_df$end)
    )
    # Apply initial chromosome style based on input parameter
    tryCatch({ seqlevelsStyle(peaks_gr) <- opt$peak_chr_style },
             error = function(e) { log_message(paste("Warning: Could not set peak_chr_style to", opt$peak_chr_style, ". Error:", e$message), "WARNING")})

}, error = function(e) {
    log_message(sprintf("Error reading or processing differential peaks CSV: %s. Details: %s", opt$diff_peaks_csv, e$message), level="ERROR")
    stop("Failed to process input CSV file.")
})
log_message(sprintf("Loaded %d differential peaks.", length(peaks_gr)))

# 2. Read GTF file and define promoter regions
log_message(sprintf("Reading GTF file from: %s", opt$gtf_file))
tryCatch({
    gtf_data <- rtracklayer::import(opt$gtf_file)
    # Apply initial chromosome style based on input parameter
    tryCatch({ seqlevelsStyle(gtf_data) <- opt$gtf_chr_style },
             error = function(e) { log_message(paste("Warning: Could not set gtf_chr_style to", opt$gtf_chr_style, ". Error:", e$message), "WARNING")})
}, error = function(e) {
    log_message(sprintf("Error reading GTF file: %s. Details: %s", opt$gtf_file, e$message), level="ERROR")
    stop("Failed to read GTF file.")
})

transcripts_gr <- gtf_data[gtf_data$type == "transcript"]
if (length(transcripts_gr) == 0) {
    log_message("No 'transcript' features found in GTF, trying 'gene' features for TSS.", level="WARNING")
    transcripts_gr <- gtf_data[gtf_data$type == "gene"]
    if (length(transcripts_gr) == 0) {
        log_message("No 'transcript' or 'gene' features found in GTF. Cannot define promoters.", level="ERROR")
        stop("Suitable features for TSS definition not found in GTF.")
    }
}
log_message(sprintf("Found %d transcript/gene features in GTF for TSS definition.", length(transcripts_gr)))

promoters_gr <- GenomicRanges::promoters(transcripts_gr, upstream=opt$promoter_upstream, downstream=opt$promoter_downstream)

if ("gene_name" %in% names(mcols(transcripts_gr))) {
    mcols(promoters_gr)$gene_name <- mcols(transcripts_gr)$gene_name
} else if ("gene_id" %in% names(mcols(transcripts_gr))) {
    mcols(promoters_gr)$gene_name <- mcols(transcripts_gr)$gene_id
} else {
    mcols(promoters_gr)$gene_name <- paste0("promoter_", seq_along(promoters_gr))
}
# Ensure gene_name is character and not factor, and handle NAs
mcols(promoters_gr)$gene_name <- as.character(mcols(promoters_gr)$gene_name)
mcols(promoters_gr)$gene_name[is.na(mcols(promoters_gr)$gene_name)] <- "UnknownGene"

log_message(sprintf("Defined %d promoter regions.", length(promoters_gr)))

# Standardize seqlevels styles for comparison (e.g., to UCSC)
# This is crucial for overlaps to work correctly if styles differ.
# We will attempt to unify them. If one is already UCSC and the other NCBI, convert NCBI to UCSC.
# If both are NCBI, convert both to UCSC for internal processing.
original_peaks_style <- seqlevelsStyle(peaks_gr)[1] # Get current style
original_promoters_style <- seqlevelsStyle(promoters_gr)[1]

common_style_to_use <- "UCSC" # Default to UCSC for internal operations
log_message(sprintf("Attempting to standardize chromosome styles. Peaks: %s, GTF/Promoters: %s. Target internal style: %s",
            original_peaks_style, original_promoters_style, common_style_to_use))

tryCatch({
    seqlevelsStyle(peaks_gr) <- common_style_to_use
}, error = function(e) {
    log_message(paste("Warning: Could not set seqlevelsStyle for peaks to", common_style_to_use, ". Error:", e$message), "WARNING")
})
tryCatch({
    seqlevelsStyle(promoters_gr) <- common_style_to_use
}, error = function(e) {
    log_message(paste("Warning: Could not set seqlevelsStyle for promoters to", common_style_to_use, ". Error:", e$message), "WARNING")
})

# Intersect seqlevels to avoid warnings/errors if one set has chromosomes not in the other
common_seqlevels <- intersect(seqlevels(peaks_gr), seqlevels(promoters_gr))
if (length(common_seqlevels) > 0) {
    peaks_gr <- keepSeqlevels(peaks_gr, common_seqlevels, pruning.mode="coarse")
    promoters_gr <- keepSeqlevels(promoters_gr, common_seqlevels, pruning.mode="coarse")
    log_message(sprintf("Kept %d common seqlevels for overlap analysis.", length(common_seqlevels)))
} else {
    log_message("No common seqlevels found between peaks and promoters after style standardization. No overlaps possible.", "WARNING")
}


# 3. Find overlaps
log_message("Finding overlaps between differential peaks and promoter regions...")
overlaps <- findOverlaps(peaks_gr, promoters_gr, ignore.strand=TRUE)

if (length(overlaps) == 0 || length(common_seqlevels) == 0) {
    log_message("No overlaps found between differential peaks and promoter regions.", level="INFO")
    file.create(opt$output_bed)
    log_message(sprintf("Empty BED file created at: %s", opt$output_bed))
    quit(save = "no", status = 0)
}
log_message(sprintf("Found %d overlaps.", length(overlaps)))

overlapping_peaks_gr <- peaks_gr[queryHits(overlaps)]
overlapping_promoters_info <- promoters_gr[subjectHits(overlaps)]
mcols(overlapping_peaks_gr)$overlapping_gene <- mcols(overlapping_promoters_info)$gene_name
mcols(overlapping_peaks_gr)$promoter_coords <- paste0(as.character(seqnames(overlapping_promoters_info)), ":", start(overlapping_promoters_info), "-", end(overlapping_promoters_info))


# 4. Format and write output BED
log_message("Formatting overlapping peaks into BED6 format.")
bed_df <- data.frame(
    chrom = as.character(seqnames(overlapping_peaks_gr)),
    chromStart = start(overlapping_peaks_gr) - 1,
    chromEnd = end(overlapping_peaks_gr),
    name = paste0(mcols(overlapping_peaks_gr)$peak_id, "_gene:", mcols(overlapping_peaks_gr)$overlapping_gene),
    score = mcols(overlapping_peaks_gr)$score,
    strand = as.character(strand(overlapping_peaks_gr))
)

bed_df$strand[bed_df$strand == "*"] <- "."

# Adjust chromosome style for output BED file
current_bed_style_is_ucsc <- all(startsWith(bed_df$chrom[!is.na(bed_df$chrom) & bed_df$chrom != ""], "chr"))

if (opt$add_chr_prefix_output) {
    if (!current_bed_style_is_ucsc) {
        log_message("Adding 'chr' prefix to chromosome names for output BED.")
        bed_df$chrom <- ifelse(!startsWith(as.character(bed_df$chrom), "chr"), paste0("chr", bed_df$chrom), as.character(bed_df$chrom))
    }
} else { # Remove "chr" prefix if not desired for output
    if (current_bed_style_is_ucsc) {
        log_message("Removing 'chr' prefix from chromosome names for output BED.")
        bed_df$chrom <- sub("^chr", "", as.character(bed_df$chrom))
    }
}

if(any(bed_df$chromStart < 0)) {
    n_negative <- sum(bed_df$chromStart < 0)
    log_message(sprintf("%d records found with negative chromStart values. Clamping to 0.", n_negative), level="WARNING")
    bed_df <- bed_df %>% dplyr::mutate(chromStart = ifelse(chromStart < 0, 0, chromStart))
}

invalid_intervals <- bed_df %>% dplyr::filter(chromStart >= chromEnd)
if (nrow(invalid_intervals) > 0) {
    log_message(sprintf("%d records found where chromStart >= chromEnd. These will be removed.", nrow(invalid_intervals)), level="WARNING")
    bed_df <- bed_df %>% dplyr::filter(chromStart < chromEnd)
}

if (nrow(bed_df) == 0) {
    log_message("No valid overlapping peaks remaining after final filtering. BED file will be empty.", level="WARNING")
    file.create(opt$output_bed)
    log_message(sprintf("Empty BED file created at: %s", opt$output_bed))
    quit(save = "no", status = 0)
}

bed_df_final <- bed_df %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand) %>%
    dplyr::distinct() # Ensure unique rows in the final output

# Sort the BED data by chromosome then start position
log_message("Sorting BED file by chromosome and start position.")
bed_df_final <- bed_df_final %>%
    dplyr::arrange(chrom, chromStart)

log_message(sprintf("Writing %d unique, sorted promoter-overlapping peaks to BED file: %s", nrow(bed_df_final), opt$output_bed))
tryCatch({
    readr::write_tsv(bed_df_final,
                     file = opt$output_bed,
                     col_names = FALSE)
}, error = function(e) {
    log_message(sprintf("Error writing output BED file: %s. Details: %s", opt$output_bed, e$message), level="ERROR")
    stop("Failed to write output BED file.")
})

log_message("Promoter overlap analysis complete.")
log_message(sprintf("Output BED file (sorted): %s", opt$output_bed))