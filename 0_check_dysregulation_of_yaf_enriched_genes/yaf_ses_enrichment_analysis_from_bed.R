#!/usr/bin/env Rscript

# YAF and SES Enrichment Analysis
# This script analyzes the enrichment of YAF vs GFP Cut&Tag data at promoters of 
# differentially expressed genes whose promoters overlap with both YAF and SES peaks.

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

required_packages <- c("DESeq2", "rtracklayer", "GenomicRanges", "ggplot2", 
                      "dplyr", "tidyr", "pheatmap", "RColorBrewer", "gridExtra",
                      "ggpubr", "scales", "viridis", "tibble", "readr") # Added readr for robust TSV reading

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE))
        BiocManager::install(pkg)
}

library(DESeq2)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(scales)
library(viridis)
library(tibble)
library(readr) # Load readr

# BASE_PATH <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
BASE_PATH <- "D:/Github/SRF_H2AK119Ub_cross_V5"

# Set paths
deseq_results_file <- file.path(BASE_PATH, "SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv")
gtf_file <- file.path(BASE_PATH, "COMMON/gencode.v43.basic.annotation.gtf") # Use relative path within project if possible

# bigwig_dir <- "F:/SRF_data/Antonella/Ub/6_bigwig"
bigwig_dir <- file.path(BASE_PATH, "DATA/6_bigwig")

# Peak file paths
yaf_peak_dir <- file.path(BASE_PATH, "DATA/5_peak_calling")
ses_peak_dir <- file.path(BASE_PATH, "SRF_SES_V5/data_from_ncbi_corrected/peaks") # Updated path for BED files

output_dir <- "YAF_SES_enrichment_results_bed" # Updated output dir name

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Function to load peaks from various formats and merge replicates
load_and_merge_peaks <- function(peak_files, format = "broadPeak") {
  all_peaks_list <- list()
  for (file in peak_files) {
    cat("  Loading peak file:", file, "\n")
    tryCatch({
      if (format == "broadPeak") {
        # Use read_delim for broadPeak to handle potential scientific notation in coordinates
        # broadPeak columns: chr, start, end, name, score, strand, signalValue, pValue, qValue
        # Read start/end as double initially to accommodate scientific notation
        col_names_broadpeak <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
        peaks_df <- read_delim(file, delim = "\t", comment = "#", 
                               col_names = col_names_broadpeak, 
                               col_types = cols(
                                 chr = col_character(),
                                 start = col_double(), # Read as double first
                                 end = col_double(),   # Read as double first
                                 name = col_character(),
                                 score = col_double(),
                                 strand = col_character(),
                                 signalValue = col_double(),
                                 pValue = col_double(),
                                 qValue = col_double()
                               ), 
                               show_col_types = FALSE) # Quiet the col spec message
        
        # Convert coordinates to integer after reading, handling potential NAs from conversion if numbers are too large for standard integer
        # Using as.integer might coerce very large numbers to NA, check if this is acceptable or if numeric is needed
        # Replace '.' strand values with '*' for GRanges compatibility
        peaks_df$strand[peaks_df$strand == "."] <- "*"
        
        # Create GRanges object first
        peaks <- GRanges(
          seqnames = peaks_df$chr,
          ranges = IRanges(start = as.integer(peaks_df$start), # Convert to integer
                           end = as.integer(peaks_df$end)),     # Convert to integer
          strand = peaks_df$strand
        )
        # Add metadata columns afterwards
        mcols(peaks) <- peaks_df[, c("name", "score", "signalValue", "pValue", "qValue")]
        # Ensure column names are set correctly
        colnames(mcols(peaks)) <- c("name", "score", "signalValue", "pValue", "qValue")

      } else if (format == "xls") {
        # Attempt to read as tab-delimited, skipping potential header lines
        # Adjust col_names and skip lines based on actual file structure
        peaks_df <- read_delim(file, delim = "\t", comment = "#", col_names = TRUE, skip = 0, show_col_types = FALSE) 
        # Assuming standard BED-like columns: chr, start, end. Adjust if needed.
        # Check column names after loading
        print(colnames(peaks_df)) 
        # Example: Adapt based on printed column names
        # Use the actual column names identified from the print statement
        # Assuming 'chr', 'start', 'end' exist based on previous run output
        required_cols <- c("chr", "start", "end") 
        if (!all(required_cols %in% colnames(peaks_df))) {
           stop("Required columns ('chr', 'start', 'end') not found in ", file)
        }
        peaks <- GRanges(
          seqnames = peaks_df$chr,
          ranges = IRanges(start = peaks_df$start, end = peaks_df$end)
          # Add strand or other metadata if available and needed
        )
      } else {
        stop("Unsupported peak format: ", format)
      }
      
      # Ensure standard chromosome naming (optional, depends on GTF)
      # seqlevelsStyle(peaks) <- "UCSC" # Example: Set to UCSC style (chr1, chr2...)
      
      all_peaks_list[[basename(file)]] <- peaks
    }, error = function(e) {
      cat("Error loading peak file:", file, "\n")
      cat("Error message:", e$message, "\n")
      # Optionally, decide whether to stop or continue if a file fails
      # stop("Failed to load peak file: ", file) # Uncomment to stop on error
    })
  }
  
  if (length(all_peaks_list) == 0) {
    stop("No peak files were loaded successfully.")
  }
  
  # Combine peaks from all replicates
  # Use reduce to merge overlapping/adjacent peaks within the combined set
  merged_peaks <- reduce(unlist(GRangesList(all_peaks_list)))
  cat(sprintf("Loaded and merged %d peaks from %d files.\n", length(merged_peaks), length(peak_files)))
  return(merged_peaks)
}


# Function to extract signal from bigWig files for a set of regions (same as original)
extract_signal_from_bigwig <- function(bigwig_file, regions) {
  tryCatch({
    cat("  Importing bigWig file:", bigwig_file, "\n")
    bw <- rtracklayer::import(bigwig_file, format = "BigWig")
    
    cat("  Chromosomes in bigWig file:", paste(unique(seqnames(bw))[1:min(5, length(unique(seqnames(bw))))], collapse=", "), "...\n")
    cat("  First few chromosomes in regions:", paste(unique(seqnames(regions))[1:min(5, length(unique(seqnames(regions))))], collapse=", "), "...\n")
    
    bw_has_chr_prefix <- any(grepl("^chr", as.character(seqnames(bw))))
    regions_has_chr_prefix <- any(grepl("^chr", as.character(seqnames(regions))))
    
    if (bw_has_chr_prefix && !regions_has_chr_prefix) {
      cat("  Adjusting regions to add 'chr' prefix...\n")
      seqlevelsStyle(regions) <- "UCSC"
    } else if (!bw_has_chr_prefix && regions_has_chr_prefix) {
      cat("  Adjusting regions to remove 'chr' prefix...\n")
      seqlevelsStyle(regions) <- "Ensembl" # Or other non-chr style
    }
    
    cat("  Finding overlaps...\n")
    overlaps <- suppressWarnings(findOverlaps(regions, bw))
    cat("  Found", length(overlaps), "overlaps\n")
    
    if (length(overlaps) == 0) {
      cat("  WARNING: No overlaps found between regions and bigWig file!\n")
      return(rep(0, length(regions))) # Return 0 signal if no overlap
    }
    
    region_signal <- tapply(
      bw$score[subjectHits(overlaps)],
      queryHits(overlaps),
      mean,
      na.rm = TRUE # Handle potential NA scores
    )
    
    all_signals <- numeric(length(regions))
    all_signals[as.numeric(names(region_signal))] <- region_signal
    all_signals[is.na(all_signals)] <- 0 # Replace NA results with 0
    
    cat("  Signal summary: min =", min(all_signals), ", max =", max(all_signals), 
        ", mean =", mean(all_signals), ", non-zero =", sum(all_signals > 0), "/", length(all_signals), "\n")
    
    return(all_signals)
  }, error = function(e) {
    cat("Error processing bigWig file:", bigwig_file, "\n")
    cat("Error message:", e$message, "\n")
    return(rep(NA, length(regions)))
  })
}

# 1. Load DESeq2 results
cat("Loading DESeq2 results...\n")
deseq_results <- read.csv(deseq_results_file, row.names = 1)
deseq_results$gene_id_no_version <- sub("\\.[0-9]+$", "", rownames(deseq_results))
cat(sprintf("Loaded %d genes from DESeq2 results.\n", nrow(deseq_results)))

# 2. Load GTF annotation
cat("\nLoading GTF annotation...\n")
gtf <- rtracklayer::import(gtf_file)
standard_chroms_no_chr <- c(as.character(1:22), "X", "Y")
standard_chroms_with_chr <- paste0("chr", standard_chroms_no_chr)
valid_chroms <- c(standard_chroms_no_chr, standard_chroms_with_chr)
gtf <- gtf[seqnames(gtf) %in% valid_chroms]
genes <- gtf[gtf$type == "gene"]
cat(sprintf("Loaded %d genes from GTF.\n", length(genes)))

# Create a mapping from gene ID to gene symbol
gene_id_to_symbol <- data.frame(
  gene_id = genes$gene_id,
  gene_name = genes$gene_name,
  stringsAsFactors = FALSE
)
gene_id_to_symbol$gene_id_no_version <- sub("\\.[0-9]+$", "", gene_id_to_symbol$gene_id)

# 3. Define promoter regions (2kb upstream, 500bp downstream of TSS)
cat("\nDefining promoter regions...\n")
genes_standard <- genes[seqnames(genes) %in% valid_chroms]
promoters <- promoters(genes_standard, upstream = 2000, downstream = 500)
promoters$gene_id_no_version <- sub("\\.[0-9]+$", "", promoters$gene_id)
cat(sprintf("Defined %d promoter regions.\n", length(promoters)))

# 4. Load and merge peak files
cat("\nLoading YAF peaks...\n")
yaf_peak_files <- list.files(yaf_peak_dir, pattern = "YAF_.*_broad_peaks_final\\.broadPeak$", full.names = TRUE)
if (length(yaf_peak_files) == 0) stop("No YAF broadPeak files found in ", yaf_peak_dir)
yaf_peaks_merged <- load_and_merge_peaks(yaf_peak_files, format = "broadPeak")

cat("\nLoading SES peaks...\n")
ses_peak_files <- list.files(ses_peak_dir, pattern = "GSM.*_broadPeak\\.bed$", full.names = TRUE) # Updated pattern for BED files
if (length(ses_peak_files) == 0) stop("No SES broadPeak BED files found in ", ses_peak_dir)
# Using format = "broadPeak" as the BED files have the 9 columns expected by this format handler
ses_peaks_merged <- load_and_merge_peaks(ses_peak_files, format = "broadPeak") # Updated format to broadPeak

# Ensure consistent chromosome naming between promoters and peaks
# Pick a style, e.g., UCSC ("chr1") or Ensembl ("1")
# If GTF uses "1", "2", etc. and peaks use "chr1", "chr2", etc.
# seqlevelsStyle(promoters) <- "UCSC" # Or match the peak style
# Or adjust peaks: seqlevelsStyle(yaf_peaks_merged) <- seqlevelsStyle(promoters)[1]
# Check styles:
cat("Promoter seqlevels style:", seqlevelsStyle(promoters)[1], "\n")
cat("YAF peak seqlevels style:", seqlevelsStyle(yaf_peaks_merged)[1], "\n")
cat("SES peak seqlevels style:", seqlevelsStyle(ses_peaks_merged)[1], "\n")

# Make styles consistent (example: set all to UCSC)
# Adjust this based on the actual styles observed
target_style <- "UCSC" # Or choose based on majority or GTF style
seqlevelsStyle(promoters) <- target_style
seqlevelsStyle(yaf_peaks_merged) <- target_style
seqlevelsStyle(ses_peaks_merged) <- target_style
cat("Applied target style:", target_style, "\n")


# 5. Find promoters overlapping with BOTH YAF and SES peaks
cat("\nFinding promoters overlapping with YAF peaks...\n")
overlaps_yaf <- findOverlaps(promoters, yaf_peaks_merged)
promoters_with_yaf <- promoters[queryHits(overlaps_yaf)]
unique_promoters_with_yaf <- unique(promoters_with_yaf) # Keep only unique promoters
cat(sprintf("Found %d unique promoters overlapping with YAF peaks.\n", length(unique_promoters_with_yaf)))

cat("Finding promoters overlapping with SES peaks...\n")
overlaps_ses <- findOverlaps(promoters, ses_peaks_merged)
promoters_with_ses <- promoters[queryHits(overlaps_ses)]
unique_promoters_with_ses <- unique(promoters_with_ses) # Keep only unique promoters
cat(sprintf("Found %d unique promoters overlapping with SES peaks.\n", length(unique_promoters_with_ses)))

cat("Finding intersection: promoters overlapping with BOTH YAF and SES peaks...\n")
# Find the common promoters by comparing their unique identifiers (e.g., gene_id_no_version)
common_gene_ids <- intersect(unique_promoters_with_yaf$gene_id_no_version, 
                             unique_promoters_with_ses$gene_id_no_version)
promoters_intersect <- promoters[promoters$gene_id_no_version %in% common_gene_ids]
# Ensure uniqueness again, although intersect should handle it
promoters_intersect <- unique(promoters_intersect) 
cat(sprintf("Found %d unique promoters overlapping with BOTH YAF and SES peaks.\n", length(promoters_intersect)))

# Get the list of gene IDs for filtering DESeq results
intersect_gene_ids <- promoters_intersect$gene_id_no_version

# 6. Filter DESeq results for genes with overlapping promoters
cat("Filtering DESeq results based on YAF/SES promoter overlap...\n")
initial_gene_count <- nrow(deseq_results)
deseq_results_filtered <- deseq_results %>%
  filter(gene_id_no_version %in% intersect_gene_ids)
cat(sprintf("Retained %d genes out of %d after filtering.\n", nrow(deseq_results_filtered), initial_gene_count))

if (nrow(deseq_results_filtered) == 0) {
  stop("No genes found in DESeq results that match the YAF/SES promoter overlap criteria. Stopping analysis.")
}

# 7. Filter for significant DEGs (padj < 0.05 and |log2FoldChange| > 1) among the filtered genes
significant_degs <- deseq_results_filtered %>%
  filter(!is.na(padj) & padj < 0.05, abs(log2FoldChange) > 1)

# Separate into up and down-regulated genes
upregulated <- significant_degs %>% filter(log2FoldChange > 0)
downregulated <- significant_degs %>% filter(log2FoldChange < 0)

cat(sprintf("Found %d significant DEGs (padj < 0.05, |log2FC| > 1) among YAF/SES overlapping genes\n", nrow(significant_degs)))
cat(sprintf("  - %d upregulated genes (YAF > GFP)\n", nrow(upregulated)))
cat(sprintf("  - %d downregulated genes (YAF < GFP)\n", nrow(downregulated)))

if (nrow(significant_degs) == 0) {
  warning("No significant DEGs found among the YAF/SES overlapping genes. Plots might be empty or analysis limited.")
  # Decide whether to stop or continue with non-significant genes if needed
  # For now, we continue but plots might be less informative.
}

# 8. Add gene symbols to DEG results
significant_degs <- merge(significant_degs, gene_id_to_symbol, 
                         by.x = "gene_id_no_version", by.y = "gene_id_no_version", all.x = TRUE)

# 9. Get promoter GRanges for the significant DEGs identified
cat("\nCreating GRanges for significant DEG promoters...\n")
# Filter the intersection promoters list to keep only those corresponding to significant DEGs
sig_promoters_intersect <- promoters_intersect[promoters_intersect$gene_id_no_version %in% significant_degs$gene_id_no_version]

# Merge DEG info with these promoters
sig_promoters_df <- as.data.frame(sig_promoters_intersect)
sig_degs_with_promoters <- merge(significant_degs, sig_promoters_df, 
                                by.x = "gene_id_no_version", by.y = "gene_id_no_version",
                                all.x = FALSE, all.y = FALSE) # Inner join

# Split into up and down regulated
up_degs_with_promoters <- sig_degs_with_promoters %>% filter(log2FoldChange > 0)
down_degs_with_promoters <- sig_degs_with_promoters %>% filter(log2FoldChange < 0)

cat(sprintf("Mapped %d significant DEGs to their promoter regions.\n", nrow(sig_degs_with_promoters)))
cat(sprintf("  - %d upregulated genes with promoters\n", nrow(up_degs_with_promoters)))
cat(sprintf("  - %d downregulated genes with promoters\n", nrow(down_degs_with_promoters)))

# Create GRanges objects for signal extraction
if (nrow(up_degs_with_promoters) > 0) {
  up_promoters_gr <- GRanges(
    seqnames = up_degs_with_promoters$seqnames,
    ranges = IRanges(start = up_degs_with_promoters$start, end = up_degs_with_promoters$end),
    strand = up_degs_with_promoters$strand,
    gene_id = up_degs_with_promoters$gene_id_no_version,
    gene_name = up_degs_with_promoters$gene_name.x, # Adjust if names conflict
    log2FC = up_degs_with_promoters$log2FoldChange
  )
} else {
  up_promoters_gr <- GRanges() # Empty GRanges if no upregulated genes
}

if (nrow(down_degs_with_promoters) > 0) {
  down_promoters_gr <- GRanges(
    seqnames = down_degs_with_promoters$seqnames,
    ranges = IRanges(start = down_degs_with_promoters$start, end = down_degs_with_promoters$end),
    strand = down_degs_with_promoters$strand,
    gene_id = down_degs_with_promoters$gene_id_no_version,
    gene_name = down_degs_with_promoters$gene_name.x, # Adjust if names conflict
    log2FC = down_degs_with_promoters$log2FoldChange
  )
} else {
  down_promoters_gr <- GRanges() # Empty GRanges if no downregulated genes
}

# 10. Extract signal from bigWig files
cat("\nExtracting signal from bigWig files for YAF/SES DEGs...\n")

# Check if bigWig files exist
yaf_files <- file.path(bigwig_dir, paste0("YAF_", 1:3, ".bw"))
gfp_files <- file.path(bigwig_dir, paste0("GFP_", 1:3, ".bw"))

all_files_exist <- all(file.exists(c(yaf_files, gfp_files)))
if (!all_files_exist) {
  missing_files <- c(yaf_files, gfp_files)[!file.exists(c(yaf_files, gfp_files))]
  cat("WARNING: Missing bigWig files:", paste(missing_files, collapse = ", "), "\n")
  stop("Cannot proceed without all bigWig files.")
}

# Initialize empty data frames
up_signal <- data.frame()
down_signal <- data.frame()

# Process upregulated genes if any exist
if (length(up_promoters_gr) > 0) {
  cat("Processing upregulated genes...\n")
  up_signal <- data.frame(
    gene_id = up_promoters_gr$gene_id,
    gene_name = up_promoters_gr$gene_name,
    log2FC_expr = up_promoters_gr$log2FC,
    regulation = "Upregulated"
  )
  for (i in 1:length(yaf_files)) {
    cat("  Processing", basename(yaf_files[i]), "...\n")
    up_signal[[paste0("YAF_", i)]] <- extract_signal_from_bigwig(yaf_files[i], up_promoters_gr)
  }
  for (i in 1:length(gfp_files)) {
    cat("  Processing", basename(gfp_files[i]), "...\n")
    up_signal[[paste0("GFP_", i)]] <- extract_signal_from_bigwig(gfp_files[i], up_promoters_gr)
  }
  # Calculate averages and log2FC binding
  up_signal$YAF_avg <- rowMeans(up_signal[, paste0("YAF_", 1:3)], na.rm = TRUE)
  up_signal$GFP_avg <- rowMeans(up_signal[, paste0("GFP_", 1:3)], na.rm = TRUE)
  up_signal$log2FC_binding <- log2((up_signal$YAF_avg + 0.1) / (up_signal$GFP_avg + 0.1)) # Add pseudocount
} else {
  cat("No upregulated significant DEGs with YAF/SES promoter overlap found. Skipping signal extraction for upregulated.\n")
}

# Process downregulated genes if any exist
if (length(down_promoters_gr) > 0) {
  cat("Processing downregulated genes...\n")
  down_signal <- data.frame(
    gene_id = down_promoters_gr$gene_id,
    gene_name = down_promoters_gr$gene_name,
    log2FC_expr = down_promoters_gr$log2FC,
    regulation = "Downregulated"
  )
  for (i in 1:length(yaf_files)) {
    cat("  Processing", basename(yaf_files[i]), "...\n")
    down_signal[[paste0("YAF_", i)]] <- extract_signal_from_bigwig(yaf_files[i], down_promoters_gr)
  }
  for (i in 1:length(gfp_files)) {
    cat("  Processing", basename(gfp_files[i]), "...\n")
    down_signal[[paste0("GFP_", i)]] <- extract_signal_from_bigwig(gfp_files[i], down_promoters_gr)
  }
  # Calculate averages and log2FC binding
  down_signal$YAF_avg <- rowMeans(down_signal[, paste0("YAF_", 1:3)], na.rm = TRUE)
  down_signal$GFP_avg <- rowMeans(down_signal[, paste0("GFP_", 1:3)], na.rm = TRUE)
  down_signal$log2FC_binding <- log2((down_signal$YAF_avg + 0.1) / (down_signal$GFP_avg + 0.1)) # Add pseudocount
} else {
  cat("No downregulated significant DEGs with YAF/SES promoter overlap found. Skipping signal extraction for downregulated.\n")
}

# Combine data if both up and down regulated genes were processed
if (nrow(up_signal) > 0 || nrow(down_signal) > 0) {
  all_signal <- rbind(up_signal, down_signal)
  
  # Handle cases where one category might be empty
  if(nrow(all_signal) > 0) {
     all_signal$regulation <- factor(all_signal$regulation, levels = c("Upregulated", "Downregulated"))
  }

  # Save processed data
  write.csv(all_signal, file.path(output_dir, "YAF_GFP_binding_at_YAF_SES_promoters.csv"), row.names = FALSE)

  # 11. Create visualizations (only if there is data)
  cat("\nCreating visualizations...\n")

  # Check if there are at least two groups to compare for some plots
  can_compare_groups <- length(unique(all_signal$regulation)) > 1

  # 1. BOXPLOT: YAF/GFP enrichment by regulation status
  p1 <- ggplot(all_signal, aes(x = regulation, y = log2FC_binding, fill = regulation)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 1, na.rm = TRUE) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1, na.rm = TRUE) + 
    theme_bw() +
    labs(title = 'YAF/GFP Enrichment at YAF/SES Promoters', # Updated title
         y = 'log2(YAF/GFP) Binding', x = '') +
    theme(legend.position = 'none', plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
    scale_fill_manual(values = c('Downregulated' = 'blue', 'Upregulated' = 'red'))
  ggsave(file.path(plots_dir, 'YAF_GFP_enrichment_boxplot.pdf'), p1, width = 6, height = 6)

  # 2. SCATTERPLOT: YAF vs GFP signal
  p2 <- ggplot(all_signal, aes(x = GFP_avg, y = YAF_avg, color = regulation)) +
    geom_point(alpha = 0.7, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, na.rm = TRUE) + 
    theme_bw() +
    labs(title = 'YAF vs GFP Signal at YAF/SES Promoters', # Updated title
         x = 'GFP Signal', y = 'YAF Signal') +
    scale_color_manual(values = c('Downregulated' = 'blue', 'Upregulated' = 'red')) +
    theme(plot.title = element_text(size = 14, face = "bold"), axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  ggsave(file.path(plots_dir, 'YAF_vs_GFP_scatterplot.pdf'), p2, width = 8, height = 6)

  # 3. BAR PLOT: Average YAF and GFP signal with error bars
  signal_summary <- all_signal %>%
    group_by(regulation) %>%
    summarize(
      YAF_mean = mean(YAF_avg, na.rm = TRUE),
      YAF_se = sd(YAF_avg, na.rm = TRUE) / sqrt(sum(!is.na(YAF_avg))),
      GFP_mean = mean(GFP_avg, na.rm = TRUE),
      GFP_se = sd(GFP_avg, na.rm = TRUE) / sqrt(sum(!is.na(GFP_avg))),
      .groups = 'drop'
    ) %>%
    pivot_longer(cols = c(YAF_mean, GFP_mean), names_to = "condition", values_to = "mean") %>%
    mutate(
      condition = sub("_mean", "", condition),
      se = ifelse(condition == "YAF", YAF_se, GFP_se)
    )

  p3 <- ggplot(signal_summary, aes(x = regulation, y = mean, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                  position = position_dodge(width = 0.8), width = 0.25, na.rm = TRUE) +
    theme_bw() +
    labs(title = 'Average Signal Intensity at YAF/SES Promoters', # Updated title
         y = 'Mean Signal', x = '', fill = 'Condition') +
    scale_fill_manual(values = c('YAF' = 'darkred', 'GFP' = 'darkblue')) +
    theme(plot.title = element_text(size = 14, face = "bold"), axis.title = element_text(size = 12),
          axis.text = element_text(size = 10), legend.position = "top")
  ggsave(file.path(plots_dir, 'average_signal_barplot.pdf'), p3, width = 8, height = 6)

  # 4. VIOLIN PLOT: Distribution of YAF/GFP ratio
  p4 <- ggplot(all_signal, aes(x = regulation, y = log2FC_binding, fill = regulation)) +
    geom_violin(trim = FALSE, alpha = 0.7, na.rm = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, na.rm = TRUE) +
    theme_bw() +
    labs(title = 'Distribution of YAF/GFP Enrichment (YAF/SES genes)', # Updated title
         y = 'log2(YAF/GFP) Binding', x = '') +
    scale_fill_manual(values = c('Downregulated' = 'blue', 'Upregulated' = 'red')) +
    theme(legend.position = 'none', plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12), axis.text = element_text(size = 10))

  # Add p-value annotation if comparison is possible
  wilcox_test_result <- NULL
  if (can_compare_groups) {
      wilcox_test_result <- tryCatch({
          test_result <- wilcox.test(log2FC_binding ~ regulation, data = all_signal)
          p_value <- test_result$p.value
          if (!is.na(p_value)) {
              p_text <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", format(p_value, digits = 3)))
              p4 <- p4 + annotate("text", x = 1.5, y = max(all_signal$log2FC_binding, na.rm = TRUE) * 0.95, 
                                  label = p_text, size = 4)
          }
          test_result # Return the result for summary
      }, error = function(e) {
          cat("Could not perform Wilcoxon test for plot annotation:", e$message, "\n")
          return(NULL)
      })
  }
  ggsave(file.path(plots_dir, 'YAF_GFP_enrichment_violin.pdf'), p4, width = 6, height = 6)

  # 5. CORRELATION PLOT: Binding vs Expression log2FC
  has_variation <- tryCatch({
      !all(is.na(all_signal$log2FC_binding)) && length(unique(na.omit(all_signal$log2FC_binding))) > 1 &&
      !all(is.na(all_signal$log2FC_expr)) && length(unique(na.omit(all_signal$log2FC_expr))) > 1 &&
      nrow(na.omit(all_signal[, c("log2FC_binding", "log2FC_expr")])) >= 3 # Need at least 3 points for cor.test
  }, error = function(e) { FALSE })

  p5 <- NULL
  cor_result <- NULL
  if (has_variation) {
      p5 <- ggplot(all_signal, aes(x = log2FC_binding, y = log2FC_expr, color = regulation)) +
          geom_point(alpha = 0.7, na.rm = TRUE) +
          geom_smooth(method = "lm", se = TRUE, color = "black", na.rm = TRUE) +
          theme_bw() +
          labs(title = 'Correlation: Binding vs Expression (YAF/SES genes)', # Updated title
               x = 'log2(YAF/GFP) Binding', y = 'log2(YAF/GFP) Expression') +
          scale_color_manual(values = c('Downregulated' = 'blue', 'Upregulated' = 'red')) +
          theme(plot.title = element_text(size = 14, face = "bold"), axis.title = element_text(size = 12),
                axis.text = element_text(size = 10))
      
      cor_result <- tryCatch({
          cor.test(~ log2FC_binding + log2FC_expr, data = all_signal, method = "pearson")
      }, error = function(e) { NULL })

      if (!is.null(cor_result) && !is.na(cor_result$estimate)) {
          cor_text <- paste0("r = ", format(cor_result$estimate, digits = 3),
                             ", p = ", format(cor_result$p.value, digits = 3))
          # Adjust annotation position if needed
          anno_x <- ifelse(min(all_signal$log2FC_binding, na.rm = TRUE) < 0, 
                           min(all_signal$log2FC_binding, na.rm = TRUE) * 0.9, 
                           min(all_signal$log2FC_binding, na.rm = TRUE) * 1.1)
          anno_y <- max(all_signal$log2FC_expr, na.rm = TRUE) * 0.9
          p5 <- p5 + annotate("text", x = anno_x, y = anno_y, label = cor_text, size = 4, hjust = 0)
      }
      ggsave(file.path(plots_dir, 'binding_expression_correlation.pdf'), p5, width = 8, height = 6)
  } else {
      cat("Skipping correlation plot due to insufficient variation or data points.\n")
  }

  # 6. DENSITY PLOT: Distribution of YAF and GFP signal
  signal_long <- all_signal %>%
    select(gene_id, regulation, YAF_avg, GFP_avg) %>%
    pivot_longer(cols = c(YAF_avg, GFP_avg), names_to = "condition", values_to = "signal") %>%
    mutate(condition = sub("_avg", "", condition))

  p6 <- ggplot(signal_long, aes(x = signal, fill = condition)) +
    geom_density(alpha = 0.5, na.rm = TRUE) +
    facet_wrap(~ regulation) +
    theme_bw() +
    labs(title = 'Distribution of YAF and GFP Signal (YAF/SES genes)', # Updated title
         x = 'Signal Intensity', y = 'Density') +
    scale_fill_manual(values = c('YAF' = 'darkred', 'GFP' = 'darkblue')) +
    theme(plot.title = element_text(size = 14, face = "bold"), axis.title = element_text(size = 12),
          axis.text = element_text(size = 10), strip.text = element_text(size = 12, face = "bold"))
  ggsave(file.path(plots_dir, 'signal_density_plot.pdf'), p6, width = 10, height = 6)

  # Create a combined figure
  plot_list <- list(p1, p4, p3, p2)
  layout_matrix <- rbind(c(1,2), c(3,4))
  if (!is.null(p5)) {
      plot_list[[5]] <- p5
      layout_matrix <- rbind(layout_matrix, c(5,5))
  }
  
  combined_plot <- grid.arrange(
      grobs = plot_list, 
      ncol = 2, 
      layout_matrix = layout_matrix,
      top = "YAF vs GFP Binding at Promoters with YAF & SES Peaks" # Updated title
  )
  ggsave(file.path(plots_dir, 'combined_plots.pdf'), combined_plot, width = 12, height = ifelse(!is.null(p5), 10, 8))

  # 12. Perform statistical tests (Wilcoxon already done for p4)
  cat("\nPerforming statistical tests...\n")
  cat('Wilcoxon rank sum test for difference in YAF/GFP enrichment between up and down-regulated genes (YAF/SES subset):\n')
  if (!is.null(wilcox_test_result) && !is.na(wilcox_test_result$p.value)) {
    print(wilcox_test_result)
  } else if (can_compare_groups) {
    cat("Could not perform Wilcoxon test (check data and variation).\n")
  } else {
    cat("Wilcoxon test skipped (only one regulation group found).\n")
  }

  # 13. Save summary statistics
  summary_stats <- all_signal %>%
    group_by(regulation) %>%
    summarize(
      mean_YAF = mean(YAF_avg, na.rm = TRUE), mean_GFP = mean(GFP_avg, na.rm = TRUE),
      median_YAF = median(YAF_avg, na.rm = TRUE), median_GFP = median(GFP_avg, na.rm = TRUE),
      mean_log2FC_binding = mean(log2FC_binding, na.rm = TRUE), median_log2FC_binding = median(log2FC_binding, na.rm = TRUE),
      sd_log2FC_binding = sd(log2FC_binding, na.rm = TRUE),
      mean_log2FC_expr = mean(log2FC_expr, na.rm = TRUE), median_log2FC_expr = median(log2FC_expr, na.rm = TRUE),
      sd_log2FC_expr = sd(log2FC_expr, na.rm = TRUE),
      n = n(), .groups = 'drop'
    )
  write.csv(summary_stats, file.path(output_dir, 'enrichment_summary_stats.csv'), row.names = FALSE)

  # 14. Create a summary report
  sink(file.path(output_dir, 'analysis_summary.txt'))
  cat('YAF Binding Analysis Summary (Filtered for Genes with YAF & SES Promoter Peaks)\n') # Updated title
  cat('================================================================================\n\n')

  cat(sprintf('Analysis restricted to genes whose promoters overlap with both YAF peaks (from %s) and SES peaks (from %s).\n', yaf_peak_dir, ses_peak_dir))
  cat(sprintf('Total promoters defined: %d\n', length(promoters)))
  cat(sprintf('Promoters overlapping YAF peaks: %d\n', length(unique_promoters_with_yaf)))
  cat(sprintf('Promoters overlapping SES peaks: %d\n', length(unique_promoters_with_ses)))
  cat(sprintf('Promoters overlapping BOTH YAF and SES peaks: %d\n', length(promoters_intersect)))
  cat(sprintf('Corresponding unique gene IDs: %d\n\n', length(intersect_gene_ids)))
  
  cat(sprintf('Genes from intersection found in DESeq results: %d\n', nrow(deseq_results_filtered)))
  cat(sprintf('Significant DEGs (padj<0.05, |log2FC|>1) among these: %d\n\n', nrow(significant_degs)))

  cat(sprintf('Number of upregulated genes analyzed: %d\n', nrow(up_signal)))
  if (nrow(up_signal) > 0) {
    cat(sprintf('Mean YAF/GFP enrichment (log2FC) for upregulated genes: %.3f ± %.3f (mean ± SD)\n\n', 
              mean(up_signal$log2FC_binding, na.rm = TRUE), sd(up_signal$log2FC_binding, na.rm = TRUE)))
  } else { cat("No upregulated genes to analyze.\n\n") }

  cat(sprintf('Number of downregulated genes analyzed: %d\n', nrow(down_signal)))
   if (nrow(down_signal) > 0) {
    cat(sprintf('Mean YAF/GFP enrichment (log2FC) for downregulated genes: %.3f ± %.3f (mean ± SD)\n\n', 
              mean(down_signal$log2FC_binding, na.rm = TRUE), sd(down_signal$log2FC_binding, na.rm = TRUE)))
  } else { cat("No downregulated genes to analyze.\n\n") }

  cat('Wilcoxon rank sum test for difference in YAF/GFP enrichment:\n')
  if (!is.null(wilcox_test_result) && !is.na(wilcox_test_result$p.value)) {
    print(wilcox_test_result)
  } else if (can_compare_groups) {
    cat("Could not perform Wilcoxon test.\n")
  } else {
     cat("Wilcoxon test skipped (only one regulation group found).\n")
  }

  cat('\nConclusion (YAF/SES subset):\n')
  if (!is.null(wilcox_test_result) && !is.na(wilcox_test_result$p.value) && wilcox_test_result$p.value < 0.05) {
    mean_up <- mean(up_signal$log2FC_binding, na.rm = TRUE)
    mean_down <- mean(down_signal$log2FC_binding, na.rm = TRUE)
    if (!is.na(mean_up) && !is.na(mean_down) && mean_up > mean_down) {
      cat('Within the YAF/SES co-bound gene set, YAF binding is significantly higher at promoters of upregulated genes compared to downregulated genes.\n')
      cat('This suggests that for these co-target genes, YAF may function as a transcriptional activator.\n')
    } else {
      cat('Within the YAF/SES co-bound gene set, YAF binding is significantly higher at promoters of downregulated genes compared to upregulated genes.\n')
      cat('This suggests that for these co-target genes, YAF binding to promoters may be associated with gene downregulation.\n')
    }
  } else if (!can_compare_groups) {
      cat('Only one group (upregulated or downregulated) was found, so comparison is not possible.\n')
  } else {
    cat('Within the YAF/SES co-bound gene set, no significant difference in YAF binding was found between promoters of up and down-regulated genes.\n')
    cat('This suggests that for these co-target genes, YAF binding alone may not determine expression direction.\n')
  }

  cat('\nCorrelation between YAF binding and gene expression (YAF/SES subset):\n')
  if (!is.null(cor_result) && !is.na(cor_result$p.value)) {
    print(cor_result)
    if (cor_result$p.value < 0.05) {
      cat(paste0('\nThere is a significant ', ifelse(cor_result$estimate > 0, 'positive', 'negative'), 
                 ' correlation between YAF binding and gene expression (r=', 
                 format(cor_result$estimate, digits=3), ').\n'))
    } else {
      cat('\nNo significant correlation found between YAF binding and gene expression.\n')
    }
  } else {
    cat("Correlation test could not be performed or was not significant.\n")
  }
  sink()

  cat("\nAnalysis complete! Results are in", output_dir, "\n")
  cat("Visualizations are in", plots_dir, "\n")

} else {
  cat("\nNo significant DEGs found in either up or down regulated categories for the YAF/SES intersection.\n")
  cat("Analysis stopped after filtering. No plots or summary statistics generated.\n")
  # Create an empty summary file indicating no results
  sink(file.path(output_dir, 'analysis_summary.txt'))
  cat('YAF Binding Analysis Summary (Filtered for Genes with YAF & SES Promoter Peaks)\n')
  cat('================================================================================\n\n')
  cat('Analysis stopped: No significant differentially expressed genes (padj<0.05, |log2FC|>1) were found whose promoters overlap with both YAF and SES peaks.\n')
  cat(sprintf('Total promoters overlapping BOTH YAF and SES peaks: %d\n', length(promoters_intersect)))
  cat(sprintf('Corresponding unique gene IDs: %d\n', length(intersect_gene_ids)))
  cat(sprintf('Genes from intersection found in DESeq results: %d\n', nrow(deseq_results_filtered)))
  sink()
  cat("Created empty summary file:", file.path(output_dir, 'analysis_summary.txt'), "\n")
}
