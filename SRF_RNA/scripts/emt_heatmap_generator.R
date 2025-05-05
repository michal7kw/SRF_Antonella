#!/usr/bin/env Rscript

# This script generates a heatmap for a specific list of genes (e.g., EMT-associated)
# using variance-stabilized expression data from a DESeq2 analysis.

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
  library(org.Hs.eg.db)  # For human gene annotation
  library(org.Mm.eg.db)  # For mouse gene annotation
})

# --- EMT Gene List ---
# Provided list of EMT-associated genes
emt_gene_list <- unique(c(
    "CDH1", "CDH2", "VIM", "ZEB1", "SNAI1", "SNAI2", "TWIST1", "TGFB1", # Corrected TGF-beta1
    "SMAD3", "ILK", "HIF1A", "EGFR", "CTNNB1", "FOXC2", "CCN2", "EPCAM",
    "CLDN4", "SERPINE1", "ESRP1", "MMP2", "MMP9", "FZD7", "LRP6", "AKT1",
    "RHOA", "RAC1", "GSK3B", "VEGFA", "NFKB1" # Corrected VEGF to VEGFA assuming common symbol
))
# Note: Duplicates like ILK, HIF1A, EGFR, ESRP1 were removed by unique().
# Note: Corrected TGF-beta1 to TGFB1 and VEGF to VEGFA (common gene symbols). Adjust if needed.

# --- Option Parsing ---
option_list <- list(
  make_option(c("-c", "--counts"), type="character", default="../results/counts/all_samples_counts.txt",
              help="Path to the featureCounts output file (e.g., all_samples_counts.txt) [default: %default]"),
  make_option(c("-m", "--metadata"), type="character", default="../metadata.csv",
              help="Path to the metadata file (CSV format, samples as rows, 'sample' and 'condition' columns required) [default: %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="results/emt_heatmap",
              help="Output directory for the heatmap [default: %default]"),
  make_option(c("--species"), type="character", default="human",
              help="Species for gene annotation ('human' or 'mouse') [default: %default]"),
  make_option(c("--heatmap_width"), type="numeric", default=10,
              help="Width of the heatmap in inches [default: %default]"),
  make_option(c("--heatmap_height"), type="numeric", default=12,
              help="Height of the heatmap in inches [default: %default]"),
  make_option(c("--heatmap_title"), type="character", default="EMT Gene Expression Heatmap",
              help="Title for the heatmap [default: %default]"),
  make_option(c("--color_palette"), type="character", default="Set1",
              help="RColorBrewer palette for condition colors [default: %default]"),
  make_option(c("--verbose"), type="logical", default=TRUE,
              help="Print verbose output [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# --- Helper Functions ---

# Function to print messages if verbose is TRUE
log_message <- function(message) {
  if (opt$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  }
}

# Function to convert gene IDs to gene symbols (adapted from original script)
get_gene_symbols <- function(gene_ids, species = opt$species) {
  if (length(gene_ids) == 0) return(character(0))
  log_message(paste("Attempting to map", length(gene_ids), "gene IDs to symbols for species:", species))

  # Select annotation db
  annotation_db <- switch(tolower(species),
                          "human" = org.Hs.eg.db,
                          "mouse" = org.Mm.eg.db,
                          { warning("Unsupported species. Choose 'human' or 'mouse'."); return(gene_ids) })

  # Detect ID type (simplified detection)
  first_id <- na.omit(gene_ids)[1]
  key_type <- "SYMBOL" # Default guess
  if (grepl("^ENS[A-Z]*[0-9]+", first_id)) {
      key_type <- "ENSEMBL"
      # Strip version if present
      if (grepl("\\.", first_id)) {
          log_message("Detected Ensembl IDs with version. Stripping version for mapping.")
          gene_ids_map <- sub("\\.[0-9]+$", "", gene_ids)
      } else {
          gene_ids_map <- gene_ids
      }
  } else if (grepl("^[0-9]+$", first_id)) {
      key_type <- "ENTREZID"
      gene_ids_map <- gene_ids
  } else {
      # Assume SYMBOL or other ID, try mapping directly
      gene_ids_map <- gene_ids
  }
  log_message(paste("Using keytype:", key_type))

  # Map IDs
  symbols <- tryCatch({
    mapIds(annotation_db,
           keys = gene_ids_map,
           column = "SYMBOL",
           keytype = key_type,
           multiVals = "first")
  }, error = function(e) {
    log_message(paste("Mapping error with keytype", key_type, ":", e$message))
    # Try SYMBOL as a fallback if the initial attempt failed
    if (key_type != "SYMBOL") {
        log_message("Trying fallback mapping with keytype: SYMBOL")
        tryCatch({
             mapIds(annotation_db, keys = gene_ids, column = "SYMBOL", keytype = "SYMBOL", multiVals = "first")
        }, error = function(e2) {
            log_message(paste("Fallback mapping error:", e2$message))
            NULL
        })
    } else {
        NULL
    }
  })

  # Handle results
  if (is.null(symbols)) {
    log_message("Gene ID to Symbol mapping failed. Using original IDs.")
    return(gene_ids)
  }

  # Fill NAs with original IDs and ensure order
  final_symbols <- symbols
  na_indices <- is.na(final_symbols)
  final_symbols[na_indices] <- gene_ids_map[na_indices] # Use the potentially version-stripped ID if mapping failed
  names(final_symbols) <- gene_ids # Name with original IDs
  final_symbols <- final_symbols[gene_ids] # Ensure original order

  log_message(paste("Mapped", sum(!na_indices), "out of", length(gene_ids), "IDs."))
  return(final_symbols)
}


# --- Input Validation ---
if (is.null(opt$counts)) {
  stop("--counts file path must be provided.", call. = FALSE)
}
if (!file.exists(opt$counts)) {
  stop(paste("Counts file not found:", opt$counts), call. = FALSE)
}
if (is.null(opt$metadata)) {
  stop("--metadata file path must be provided.", call. = FALSE)
}
if (!file.exists(opt$metadata)) {
  stop(paste("Metadata file not found:", opt$metadata), call. = FALSE)
}
if (!tolower(opt$species) %in% c("human", "mouse")) {
  stop("Invalid species specified. Choose 'human' or 'mouse'.", call. = FALSE)
}

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", opt$output_dir))
}

# --- Load Data ---

# Read full count matrix
log_message(paste("Reading full count matrix from:", opt$counts))
raw_counts <- read.delim(opt$counts, comment.char = "#", header = TRUE, sep = "\t", check.names = FALSE)

# Set Geneid as row names and select only count columns (7th column onwards)
# Assumes standard featureCounts output format
if (ncol(raw_counts) < 7 || !("Geneid" %in% colnames(raw_counts))) {
    stop("Counts file does not seem to be in standard featureCounts format (missing 'Geneid' or < 7 columns).")
}
rownames(raw_counts) <- raw_counts$Geneid
count_matrix_full <- raw_counts[, 7:ncol(raw_counts), drop = FALSE]

# Clean up sample names in the count matrix (remove path prefixes if present)
# Example: results/star/C1/C1_Aligned.sortedByCoord.out.bam -> C1
# Adjust regex if your sample names have a different structure
colnames(count_matrix_full) <- sub(".*/(.*?)/\\1_Aligned.*", "\\1", colnames(count_matrix_full))
colnames(count_matrix_full) <- sub("_Aligned.*", "", colnames(count_matrix_full)) # More general cleanup
colnames(count_matrix_full) <- sub("\\.[sb]am$", "", colnames(count_matrix_full)) # Remove .bam/.sam suffix
log_message(paste("Cleaned sample names:", paste(head(colnames(count_matrix_full)), collapse=", ")))

# Read metadata
log_message(paste("Reading metadata from:", opt$metadata))
metadata_full <- read.csv(opt$metadata, row.names = "sample", check.names = FALSE)
log_message(paste("Metadata loaded for", nrow(metadata_full), "samples."))

# Ensure samples in count matrix and metadata are aligned
samples_in_counts <- colnames(count_matrix_full)
samples_in_metadata <- rownames(metadata_full)
common_samples <- intersect(samples_in_counts, samples_in_metadata)

if(length(common_samples) == 0) {
  stop("No common samples found between metadata and count matrix. Check sample names.")
}
log_message(paste("Found", length(common_samples), "common samples between counts and metadata."))

if(length(common_samples) < length(samples_in_counts)) {
  log_message("Warning: Subsetting count matrix to samples present in metadata.")
  count_matrix_full <- count_matrix_full[, common_samples, drop = FALSE]
}
if(length(common_samples) < length(samples_in_metadata)) {
  log_message("Warning: Subsetting metadata to samples present in count matrix.")
  metadata_full <- metadata_full[common_samples, , drop = FALSE]
}

# Order metadata rows to match count matrix columns
metadata_full <- metadata_full[colnames(count_matrix_full), , drop = FALSE]

# Verify that sample order matches before filtering
if (!all(rownames(metadata_full) == colnames(count_matrix_full))) {
  stop("Initial sample order mismatch after filtering/ordering metadata and counts. This should not happen.")
}

# --- Filter for YAF and GFP conditions ---
target_conditions_lower <- tolower(c("YAF", "GFP")) # Use lower case for comparison
log_message(paste("Filtering samples for conditions (case-insensitive):", paste(c("YAF", "GFP"), collapse=", ")))

# Ensure 'condition' column exists
if (!"condition" %in% colnames(metadata_full)) {
    stop("Metadata file must contain a 'condition' column.")
}

# Perform case-insensitive filtering
keep_rows <- tolower(trimws(metadata_full$condition)) %in% target_conditions_lower
samples_to_keep <- rownames(metadata_full)[keep_rows]

if (length(samples_to_keep) == 0) {
    stop(paste("No samples found matching conditions (case-insensitive):", paste(c("YAF", "GFP"), collapse=", ")))
}
log_message(paste("Keeping", length(samples_to_keep), "samples for analysis:", paste(samples_to_keep, collapse=", ")))

# Filter metadata and count matrix using the identified samples
log_message(paste("Subsetting count matrix and metadata for", length(samples_to_keep), "samples."))
metadata_filtered <- metadata_full[samples_to_keep, , drop = FALSE]
count_matrix_filtered <- count_matrix_full[, samples_to_keep, drop = FALSE]

# Log the conditions actually present after filtering
final_conditions_in_metadata <- unique(metadata_filtered$condition)
log_message(paste("Conditions present in filtered metadata:", paste(final_conditions_in_metadata, collapse=", ")))
log_message(paste("Samples present in filtered count matrix:", paste(colnames(count_matrix_filtered), collapse=", ")))


# Verify order again after filtering
if (!all(rownames(metadata_filtered) == colnames(count_matrix_filtered))) {
  stop("Sample order mismatch after filtering for target conditions. This should not happen.")
}
# --- End Filtering ---


# Create DESeqDataSet using only filtered data
log_message("Creating DESeqDataSet using filtered YAF/GFP samples")
dds_filtered <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_filtered,
  design = ~ condition # Design based on the conditions present
)

# Filter out low-count genes before VST (optional but recommended)
dds_filtered <- dds_filtered[rowSums(counts(dds_filtered)) >= 10, ]
log_message(paste("Filtered low-count genes. Kept", nrow(dds_filtered), "genes."))

# Get variance stabilized transformed data using only filtered samples
log_message("Performing variance stabilizing transformation on filtered YAF/GFP samples")
vst_data <- vst(dds_filtered, blind = FALSE)

log_message(paste("VST computed for", ncol(vst_data), "samples and", nrow(vst_data), "genes."))
# --- Verify Samples in VST Data ---
vst_samples <- colnames(vst_data)
log_message(paste("Samples present in the final VST object:", paste(vst_samples, collapse=", ")))
# Cross-reference with filtered metadata
if (!all(vst_samples %in% rownames(metadata_filtered))) {
    warning("Mismatch between VST samples and filtered metadata rownames!")
}
vst_conditions <- colData(vst_data)$condition
log_message(paste("Conditions present in the final VST object colData:", paste(unique(vst_conditions), collapse=", ")))
if (!all(tolower(unique(vst_conditions)) %in% tolower(c("YAF", "GFP")))) {
     warning("Unexpected conditions found in VST colData after filtering!")
}
# --- End Verification ---


# --- Prepare Data for Heatmap ---

vst_assay <- assay(vst_data)
vst_gene_ids <- rownames(vst_assay)

# Get gene symbols for the VST data rows
log_message("Mapping VST gene IDs to symbols.")
vst_gene_symbols <- get_gene_symbols(vst_gene_ids, species = opt$species)
names(vst_gene_symbols) <- vst_gene_ids # Keep original IDs as names

# Find intersection between EMT list and genes present in VST data (using symbols)
emt_genes_in_data_symbols <- intersect(emt_gene_list, vst_gene_symbols)
log_message(paste("Found", length(emt_genes_in_data_symbols), "out of", length(emt_gene_list), "EMT genes in the VST data."))

if (length(emt_genes_in_data_symbols) < 2) {
  stop("Fewer than 2 EMT genes found in the expression data. Cannot generate heatmap.", call. = FALSE)
}

# Get the original IDs corresponding to the found symbols
emt_genes_in_data_ids <- names(vst_gene_symbols)[match(emt_genes_in_data_symbols, vst_gene_symbols)]
# Handle potential multiple IDs mapping to the same symbol (take the first match)
emt_genes_in_data_ids <- emt_genes_in_data_ids[!is.na(emt_genes_in_data_ids)] # Remove NAs if any symbol wasn't found back in names
emt_genes_in_data_symbols <- vst_gene_symbols[emt_genes_in_data_ids] # Get the symbols corresponding to the selected IDs

# Filter the VST assay matrix
emt_assay_subset <- vst_assay[emt_genes_in_data_ids, , drop = FALSE]
rownames(emt_assay_subset) <- emt_genes_in_data_symbols # Use symbols for heatmap rows

# Scale the data
log_message("Scaling expression data for heatmap.")
scaled_data <- t(scale(t(emt_assay_subset)))
scaled_data[is.nan(scaled_data)] <- 0 # Handle potential NaNs if variance is zero

# Prepare annotation using the *filtered* metadata
log_message("Preparing heatmap annotation using filtered metadata.")
annotation_col <- data.frame(Condition = metadata_filtered$condition, # Use the filtered metadata
                             row.names = rownames(metadata_filtered)) # Use rownames from filtered metadata

# Verify conditions in the annotation dataframe
final_conditions <- unique(annotation_col$Condition)
log_message(paste("Conditions included in final heatmap annotation:", paste(final_conditions, collapse=", ")))
if (!all(final_conditions %in% c("YAF", "GFP"))) {
    warning("Unexpected conditions found in the final annotation data. Check filtering logic.")
    # Filter annotation_col again just in case
    annotation_col <- annotation_col[annotation_col$Condition %in% c("YAF", "GFP"), , drop = FALSE]
    final_conditions <- unique(annotation_col$Condition)
    log_message(paste("Conditions after re-filtering annotation:", paste(final_conditions, collapse=", ")))
}


# Prepare annotation colors based *only* on the conditions present in the filtered data
conditions <- final_conditions # Use the verified conditions
num_colors_needed <- length(conditions)

# Get the full palette
full_palette <- brewer.pal(brewer.pal.info[opt$color_palette, "maxcolors"], opt$color_palette)

# Select only the required number of colors from the start of the palette
if (num_colors_needed > length(full_palette)) {
    warning(paste("Number of conditions (", num_colors_needed, ") exceeds colors available in palette", opt$color_palette, ". Colors will be recycled."))
    selected_colors <- rep(full_palette, length.out = num_colors_needed)
} else {
    selected_colors <- full_palette[1:num_colors_needed]
}

# Assign colors to conditions
condition_colors <- setNames(selected_colors, conditions)
log_message("Assigned colors to conditions:")
log_message(paste(names(condition_colors), condition_colors, sep = "=", collapse = ", "))
ann_colors <- list(Condition = condition_colors)

# --- Generate Heatmap ---

log_message("Generating heatmap.")

# Define output files
heatmap_file_pdf <- file.path(opt$output_dir, "emt_gene_heatmap.pdf")
heatmap_file_png <- file.path(opt$output_dir, "emt_gene_heatmap.png")

# Define color scale and breaks (dynamic based on scaled data)
lim <- max(abs(scaled_data), na.rm = TRUE)
lim <- ceiling(lim * 10) / 10 # Round up slightly
if (lim == 0) lim <- 1 # Avoid zero range
my_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
my_breaks <- seq(-lim, lim, length.out = 101)

# Create heatmap function call (to avoid repetition)
create_heatmap <- function() {
    pheatmap(scaled_data,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = my_colors,
             breaks = my_breaks,
             border_color = NA, # "grey60" or NA
             main = opt$heatmap_title,
             clustering_method = "ward.D2",
             show_colnames = TRUE, # Show sample names
             show_rownames = TRUE, # Show gene symbols
             fontsize = 10,
             fontsize_row = max(6, 10 - nrow(scaled_data) / 10), # Adjust row font size based on number of genes
             fontsize_col = 8)
}

# Save PDF
tryCatch({
    pdf(heatmap_file_pdf, width = opt$heatmap_width, height = opt$heatmap_height)
    create_heatmap()
    dev.off()
    log_message(paste("Heatmap saved to:", heatmap_file_pdf))
}, error = function(e) {
    warning(paste("Failed to save PDF heatmap:", e$message))
    if (exists("dev.off")) dev.off() # Close device if error occurred during plotting
})


# Save PNG
tryCatch({
    png(heatmap_file_png, width = opt$heatmap_width, height = opt$heatmap_height, units = "in", res = 300)
    create_heatmap()
    dev.off()
    log_message(paste("Heatmap saved to:", heatmap_file_png))
}, error = function(e) {
    warning(paste("Failed to save PNG heatmap:", e$message))
     if (exists("dev.off")) dev.off() # Close device if error occurred during plotting
})


log_message("Script finished.")