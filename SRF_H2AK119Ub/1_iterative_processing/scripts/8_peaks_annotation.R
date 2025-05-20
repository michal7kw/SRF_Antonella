# Filename: 8_annotation_and_piechart_per_sample.R
#
# Description:
# This script performs peak annotation for individual BED files (e.g., consensus peaks
# for different samples) and generates a detailed pie chart of genomic feature
# distribution for each sample.
#
# Input files (hardcoded for this example):
# - ./analysis/6_consensus_peaks/GFP_consensus_peaks.bed
# - ./analysis/6_consensus_peaks/YAF_consensus_peaks.bed
#
# Output files (example for GFP sample, YAF will be similar):
# In ./analysis/8_annotation_and_enrichment_per_sample/GFP_annotation/:
#   figures/
#     - GFP_detailed_pie_chart.pdf: Detailed pie chart with genomic feature distribution for GFP
#   tables/
#     - GFP_peak_annotation.csv: Detailed peak annotations for GFP
#   GFP_peak_annotation.rds: R object with full annotation data for GFP
#
# Dependencies:
# - ChIPseeker for peak annotation
# - org.Hs.eg.db for gene ID mapping (though not directly used for pie chart, good for consistency)
# - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
# - tidyverse for data manipulation
# - ggplot2 for custom visualizations
# - GenomicRanges for GRanges objects

# Load necessary libraries
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(org.Hs.eg.db) # Kept for potential extensions, though not strictly needed for this pie chart
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(tidyverse)
    library(GenomicRanges)
    library(ggplot2)
    library(rtracklayer) # For importing BED files
})

# Attempt to source utility functions if available (modify path if needed)
if (file.exists("scripts/utils.R")) {
    source("scripts/utils.R")
} else {
    # Define dummy functions if utils.R is not found, to avoid errors
    log_message <- function(message, level = "INFO") {
        cat(sprintf("[%s %s] %s\n", Sys.time(), level, message))
    }
    create_dirs <- function(dirs) {
        for (d in dirs) {
            if (!dir.exists(d)) {
                dir.create(d, recursive = TRUE)
                log_message(sprintf("Created directory: %s", d))
            }
        }
    }
}

#' Create a detailed pie chart showing genomic feature distribution
#' @param anno_data Data frame with annotation data (must contain 'annotation' column)
#' @param title Optional title for the pie chart
#' @return ggplot object with detailed pie chart
create_detailed_pie_chart <- function(anno_data, title = NULL) {
    # Make a copy of the annotation data
    anno_df <- anno_data # In this version, anno_data is directly from as.data.frame(peak_anno)
    
    # Create a more detailed categorization
    anno_df$detailed_annotation <- NA
    
    # Promoter regions with different distances
    anno_df$detailed_annotation[grepl("Promoter \\(<=1kb\\)", anno_df$annotation)] <- "Promoter (<=1kb)"
    anno_df$detailed_annotation[grepl("Promoter \\(1-2kb\\)", anno_df$annotation)] <- "Promoter (1-2kb)"
    anno_df$detailed_annotation[grepl("Promoter \\(2-3kb\\)", anno_df$annotation)] <- "Promoter (2-3kb)"
    
    # UTR regions
    anno_df$detailed_annotation[grepl("5' UTR", anno_df$annotation)] <- "5' UTR"
    anno_df$detailed_annotation[grepl("3' UTR", anno_df$annotation)] <- "3' UTR"
    
    # Exon and Intron regions
    anno_df$detailed_annotation[grepl("Exon", anno_df$annotation)] <- "Exon" # Handles "Exon (coding)" etc.
    anno_df$detailed_annotation[grepl("Intron", anno_df$annotation)] <- "Intron" # Handles "Intron (NR_...)" etc.
    
    # Downstream and intergenic
    anno_df$detailed_annotation[grepl("Downstream", anno_df$annotation)] <- "Downstream (<=300bp)" # ChIPseeker default for "Downstream" often refers to <3kb, this script previously used 300bp label
    anno_df$detailed_annotation[grepl("Distal Intergenic", anno_df$annotation)] <- "Distal Intergenic"
    
    # Handle cases where detailed_annotation might still be NA (e.g. "Other")
    # For simplicity, we group remaining NAs or less specific annotations if needed,
    # or ensure the grepl patterns are comprehensive.
    # Here, we will rely on the existing categories.
    # If specific categories like "Other" from ChIPseeker output are desired, they can be added.
    # For now, anything not matching the above will be NA and might be implicitly excluded or cause issues if not handled.
    # Let's refine to ensure all annotations get a category, or are grouped into 'Other'.
    
    # Re-categorize based on ChIPseeker's typical output structure more robustly
    anno_df <- anno_df %>%
      mutate(detailed_annotation = case_when(
        grepl("Promoter \\(<=1kb\\)", annotation) ~ "Promoter (<=1kb)",
        grepl("Promoter \\(1-2kb\\)", annotation) ~ "Promoter (1-2kb)",
        grepl("Promoter \\(2-3kb\\)", annotation) ~ "Promoter (2-3kb)",
        grepl("5' UTR", annotation) ~ "5' UTR",
        grepl("3' UTR", annotation) ~ "3' UTR",
        grepl("Exon", annotation) ~ "Exon", # Covers 1st Exon, Other Exon
        grepl("Intron", annotation) ~ "Intron", # Covers 1st Intron, Other Intron
        grepl("Downstream", annotation) ~ "Downstream (gene end)", # Default ChIPseeker label
        grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
        TRUE ~ "Other" # Catch-all for annotations not fitting above
      ))

    # Calculate percentages
    anno_summary <- anno_df %>%
        group_by(detailed_annotation) %>%
        summarise(count = n(), .groups = 'drop') %>%
        mutate(percentage = count / sum(count) * 100) %>%
        filter(!is.na(detailed_annotation)) %>% # Ensure NA annotations are not plotted
        arrange(desc(percentage))
    
    # Define colors similar to the example image
    color_palette <- c(
        "Promoter (<=1kb)" = "#A6CEE3",
        "Promoter (1-2kb)" = "#1F78B4",
        "Promoter (2-3kb)" = "#B2DF8A",
        "5' UTR" = "#33A02C",
        "3' UTR" = "#FB9A99",
        "Exon" = "#E31A1C",
        "Intron" = "#FF7F00",
        "Downstream (gene end)" = "#CAB2D6", # Adjusted from original to match ChIPseeker typical category
        "Distal Intergenic" = "#6A3D9A", # Adjusted from original for consistency
        "Other" = "#FFFF99" # Color for 'Other' category
    )
    
    # Create the pie chart
    pie_chart <- ggplot(anno_summary, aes(x = "", y = percentage, fill = detailed_annotation)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = color_palette, name = "Genomic Feature",
                          limits = names(color_palette)[names(color_palette) %in% anno_summary$detailed_annotation]) + # Ensure legend order and presence
        theme_void() +
        theme(legend.position = "right") +
        geom_text(aes(label = ifelse(percentage >= 2, paste0(round(percentage, 1), "%"), "")), # Show label if >= 2%
                  position = position_stack(vjust = 0.5),
                  size = 3, check_overlap = TRUE)
    
    # Add title if provided
    if (!is.null(title)) {
        pie_chart <- pie_chart + ggtitle(title)
    }
    
    return(pie_chart)
}

#' Standardize chromosome names to UCSC style and filter valid chromosomes
#' @param gr GRanges object to standardize
#' @return GRanges object with standardized chromosome names
standardize_chromosomes <- function(gr) {
    # Define standard chromosomes (1-22, X, Y, M if needed)
    # Keeping M (mitochondrial) might be relevant depending on the study
    standard_chroms_ucsc <- paste0("chr", c(1:22, "X", "Y", "M"))
    
    current_seqlevels <- seqlevels(gr)
    new_seqlevels <- current_seqlevels
    
    # Add 'chr' prefix if missing for numeric, X, Y, M chromosomes
    needs_prefix <- grep("^([1-9]|1[0-9]|2[0-2]|[XYM])$", current_seqlevels)
    if (length(needs_prefix) > 0) {
      new_seqlevels[needs_prefix] <- paste0("chr", current_seqlevels[needs_prefix])
    }
    
    # Rename seqlevels in the object
    tryCatch({
        gr <- renameSeqlevels(gr, new_seqlevels)
    }, error = function(e) {
        log_message(paste("Error renaming seqlevels:", e$message), level = "WARNING")
        # If renaming fails due to duplicate names (e.g. 'chr1' and '1' becoming 'chr1'),
        # it means some might already be in UCSC style. We proceed with filtering.
    })

    # Filter to keep only standard chromosomes
    gr_filtered <- gr[seqnames(gr) %in% standard_chroms_ucsc]
    
    # Update sequence levels to only those present in the filtered GRanges and are standard
    final_seqlevels <- intersect(standard_chroms_ucsc, seqlevels(gr_filtered))
    seqlevels(gr_filtered) <- final_seqlevels # Drop unused levels
    
    # Ensure genome information is consistent if available
    # If your BED files have genome info, it might be useful to check/set it here.
    # For now, assuming hg38 as per TxDb.
    
    return(gr_filtered)
}


#' Perform peak annotation and generate a detailed pie chart for a single sample
#' @param bed_file_path Path to the input BED file
#' @param sample_name A string identifier for the sample (e.g., "GFP", "YAF")
#' @param output_dir_base Base directory for output files
#' @return Invisibly returns the annotation object
annotate_and_generate_pie_for_sample <- function(bed_file_path, sample_name, output_dir_base) {
    log_message(sprintf("Starting annotation for sample: %s from file: %s", sample_name, bed_file_path))

    # Define output directories
    sample_output_dir <- file.path(output_dir_base, paste0(sample_name, "_annotation"))
    figures_dir <- file.path(sample_output_dir, "figures")
    tables_dir <- file.path(sample_output_dir, "tables")

    # Create output directories
    create_dirs(c(figures_dir, tables_dir))

    # Load peaks from BED file
    log_message(sprintf("Reading peaks from BED file: %s", bed_file_path))
    if (!file.exists(bed_file_path)) {
        log_message(sprintf("ERROR: BED file not found: %s", bed_file_path), level = "ERROR")
        return(invisible(NULL))
    }
    # peaks_gr <- ChIPseeker::readPeakFile(bed_file_path, as = "GRanges") # ChIPseeker's function
    peaks_gr <- tryCatch({
        rtracklayer::import(bed_file_path, format = "BED")
    }, error = function(e) {
        log_message(sprintf("ERROR: Failed to read BED file %s: %s", bed_file_path, e$message), level = "ERROR")
        return(NULL)
    })

    if (is.null(peaks_gr) || length(peaks_gr) == 0) {
        log_message(sprintf("No peaks loaded from %s or file is empty. Skipping.", bed_file_path), level = "WARNING")
        return(invisible(NULL))
    }
    log_message(sprintf("Loaded %d peaks for sample %s.", length(peaks_gr), sample_name))

    # Standardize chromosome names
    peaks_gr_std <- standardize_chromosomes(peaks_gr)
    if (length(peaks_gr_std) == 0) {
        log_message(sprintf("No peaks remaining after chromosome standardization for %s. Skipping.", sample_name), level = "WARNING")
        return(invisible(NULL))
    }
    log_message(sprintf("%d peaks remaining after chromosome standardization for %s.", length(peaks_gr_std), sample_name))

    # Annotate peaks
    log_message(sprintf("Performing peak annotation for %s...", sample_name))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno <- annotatePeak(peaks_gr_std,
        TxDb = txdb,
        tssRegion = c(-3000, 3000), # Define promoter region
        verbose = FALSE
    )

    if (is.null(peak_anno) || length(peak_anno@anno) == 0) {
        log_message(sprintf("Peak annotation failed or resulted in no annotations for %s.", sample_name), level = "WARNING")
        return(invisible(NULL))
    }
    
    # Save annotation results
    peak_annotation_rds_file <- file.path(sample_output_dir, paste0(sample_name, "_peak_annotation.rds"))
    saveRDS(peak_anno, peak_annotation_rds_file)

    peak_annotation_csv_file <- file.path(tables_dir, paste0(sample_name, "_peak_annotation.csv"))
    write.csv(as.data.frame(peak_anno),
        peak_annotation_csv_file,
        row.names = FALSE
    )
    log_message(sprintf("Saved peak annotation RDS and CSV for %s.", sample_name))

    # Prepare data for pie chart
    # The create_detailed_pie_chart function expects a data frame with an 'annotation' column.
    # It was originally designed to also handle 'fold_change' for differential analysis,
    # but that part is conditional and won't be used if 'direction' is NULL.
    annotation_df_for_pie <- as.data.frame(peak_anno@anno) # peak_anno is a csAnno object

    # Create and save detailed pie chart
    log_message(sprintf("Generating detailed pie chart for %s...", sample_name))
    detailed_pie <- create_detailed_pie_chart(annotation_df_for_pie, 
                                              title = paste("Genomic Feature Distribution -", sample_name))
    
    detailed_pie_pdf_path <- file.path(figures_dir, paste0(sample_name, "_detailed_pie_chart.pdf"))
    pdf(detailed_pie_pdf_path, width = 10, height = 8)
    print(detailed_pie)
    dev.off()
    log_message(sprintf("Saved detailed pie chart for %s to: %s", sample_name, detailed_pie_pdf_path))

    log_message(sprintf("Completed annotation for sample: %s", sample_name))
    return(invisible(peak_anno))
}

# --- Main Execution Block ---

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript 8_peaks_annotation.R <input_dir> <output_dir_base>", call. = FALSE)
}

INPUT_DIR_BASE <- args[1]
PIPELINE_OUTPUT_DIR_BASE <- args[2]

log_message(sprintf("Input directory set to: %s", INPUT_DIR_BASE))
log_message(sprintf("Output directory base set to: %s", PIPELINE_OUTPUT_DIR_BASE))

log_message("Starting per-sample peak annotation and pie chart generation pipeline...")

create_dirs(PIPELINE_OUTPUT_DIR_BASE)

# Define input BED files and sample names
# Construct full paths based on the input directory argument
input_files <- list(
    GFP = file.path(INPUT_DIR_BASE, "GFP_consensus_peaks.bed"),
    YAF = file.path(INPUT_DIR_BASE, "YAF_consensus_peaks.bed")
)

# Process each sample
for (sample_id in names(input_files)) {
    bed_file <- input_files[[sample_id]]
    
    if (!file.exists(bed_file)) {
        log_message(sprintf("WARNING: Input BED file not found for sample %s: %s. Skipping.", sample_id, bed_file), level = "WARNING")
        next # Skip to the next sample
    }
    
    annotate_and_generate_pie_for_sample(
        bed_file_path = bed_file,
        sample_name = sample_id,
        output_dir_base = PIPELINE_OUTPUT_DIR_BASE
    )
}

log_message("Pipeline completed successfully.")