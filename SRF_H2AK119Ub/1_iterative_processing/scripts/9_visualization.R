# This script creates visualizations for ChIP-seq peak analysis results.
# It generates plots and statistics to help interpret the differential binding analysis
# between YAF and GFP samples.
#
# Input files:
# - analysis/7_differential_binding/diffbind_broad/significant_peaks.rds: GRanges object with differential peaks
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_full.csv: Annotated gene list
#
# Output files in analysis/9_visualization/{peak_type}/:
#   peak_analysis/
#     - peak_distribution.pdf: Bar plot showing distribution of peak categories
#     - ma_plot_enhanced.pdf: MA plot with fold change vs concentration
#     - volcano_plot_enhanced.pdf: Volcano plot of significance vs fold change
#   gene_analysis/
#     - tss_distribution.pdf: Histogram of peak distances to TSS
#   summary_statistics/
#     - summary_stats.csv: Table with key analysis metrics
#
# Dependencies:
# - ggplot2 for plotting
# - gridExtra for plot arrangement
# - scales for axis formatting
# - utils.R for helper functions

# Load necessary libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(scales)
    library(ChIPseeker)
    library(GenomicRanges)
})

source("scripts/utils.R")

#' Create comprehensive visualizations for peak analysis
#' @param output_dir Base directory for output files
#' @param peak_type Character, "broad"
#' @param diff_peaks DiffBind results as GRanges or data.frame
#' @param gene_list Data frame with gene annotations including distance to TSS
#' @return List containing all generated plots and summary statistics
create_visualizations <- function(output_dir, peak_type, diff_peaks, gene_list) {
    log_message(sprintf("Creating visualizations for %s peaks", peak_type))
    
    # Define base output directory
    base_dir <- file.path(output_dir, peak_type)
    
    # Create output directories with full paths
    peak_analysis_dir <- file.path(base_dir, "peak_analysis")
    gene_analysis_dir <- file.path(base_dir, "gene_analysis")
    summary_statistics_dir <- file.path(base_dir, "summary_statistics")
    
    dirs <- c(peak_analysis_dir, gene_analysis_dir, summary_statistics_dir)
    create_dirs(dirs)
    
    # Convert peaks to data frame if needed
    if (is(diff_peaks, "GRanges")) {
        diff_peaks <- as.data.frame(diff_peaks)
    }
    
    # 1. Peak Distribution Plot - Shows counts of total, YAF-enriched and GFP-enriched peaks
    log_message("Creating peak distribution plot...")
    peak_dist <- data.frame(
        Category = c("Total", "YAF-enriched", "GFP-enriched"),
        Count = c(
            nrow(diff_peaks),
            sum(diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05),
            sum(diff_peaks$Fold < 0 & diff_peaks$FDR < 0.05)
        )
    )
    
    peak_dist_plot <- ggplot(peak_dist, aes(x = reorder(Category, -Count), y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
        theme_minimal() +
        labs(x = "Category",
             y = "Number of Peaks",
             title = paste("Distribution of", peak_type, "Peaks")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_plot(peak_dist_plot,
              file.path(peak_analysis_dir, "peak_distribution.pdf"))
    
    # 2. Enhanced MA Plot - Shows relationship between fold change and mean concentration
    log_message("Creating enhanced MA plot...")
    ma_plot <- ggplot(diff_peaks, aes(x = log2(Conc), y = Fold)) +
        geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
        scale_color_manual(values = c("grey50", "red3")) +
        scale_size_continuous(range = c(0.5, 2)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        labs(x = "log2 Mean Concentration",
             y = "log2 Fold Change (YAF/GFP)",
             title = paste("MA Plot:", peak_type, "Peaks")) +
        theme_minimal() +
        theme(legend.position = "right")
    
    save_plot(ma_plot,
              file.path(peak_analysis_dir, "ma_plot_enhanced.pdf"))
    
    # 3. Enhanced Volcano Plot - Shows significance vs fold change
    log_message("Creating enhanced volcano plot...")
    volcano_plot <- ggplot(diff_peaks, aes(x = Fold, y = -log10(FDR))) +
        geom_point(aes(color = FDR < 0.05 & abs(Fold) > 1,
                       size = abs(Fold)), alpha = 0.6) +
        scale_color_manual(values = c("grey70", "red3")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
        labs(x = "log2 Fold Change (YAF/GFP)",
             y = "-log10(FDR)",
             title = paste("Volcano Plot:", peak_type, "Peaks")) +
        theme_minimal()
    
    save_plot(volcano_plot,
              file.path(peak_analysis_dir, "volcano_plot_enhanced.pdf"))
    
    # 4. Gene Distance to TSS Distribution - Shows peak locations relative to genes
    log_message("Creating TSS distance distribution plot...")
    tss_dist_plot <- ggplot(gene_list, aes(x = distanceToTSS)) +
        geom_histogram(bins = 50, fill = "darkred", color = "black", alpha = 0.7) +
        theme_minimal() +
        labs(x = "Distance to TSS (bp)",
             y = "Frequency",
             title = paste("Distribution of Peak Distances to TSS -", peak_type)) +
        scale_x_continuous(labels = scales::comma)
    
    save_plot(tss_dist_plot,
              file.path(gene_analysis_dir, "tss_distribution.pdf"))
    
    # 5. Create summary statistics - Key metrics about the analysis
    log_message("Generating summary statistics...")
    summary_stats <- data.frame(
        Metric = c(
            "Total Peaks",
            "Significant Peaks (FDR < 0.05)",
            "YAF-enriched Peaks",
            "GFP-enriched Peaks",
            "Unique Genes",
            "Median Distance to TSS"
        ),
        Value = c(
            nrow(diff_peaks),
            sum(diff_peaks$FDR < 0.05),
            sum(diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05),
            sum(diff_peaks$Fold < 0 & diff_peaks$FDR < 0.05),
            length(unique(gene_list$SYMBOL)),
            median(abs(gene_list$distanceToTSS))
        )
    )
    
    write.csv(summary_stats,
              file.path(summary_statistics_dir, "summary_stats.csv"),
              row.names = FALSE)
    
    log_message(sprintf("Completed visualizations for %s peaks", peak_type))
    return(list(
        peak_dist_plot = peak_dist_plot,
        ma_plot = ma_plot,
        volcano_plot = volcano_plot,
        tss_dist_plot = tss_dist_plot,
        summary_stats = summary_stats
    ))
}

# Get output directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Output directory argument is required")
}
output_dir <- args[1]

# Process both peak types
log_message("Starting visualization process...")

# Process broad peaks
broad_peaks_file <- "analysis/7_differential_binding/diffbind_broad/significant_peaks.rds"
broad_genes_file <- "analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_full.csv"

log_message(sprintf("Loading broad peaks from %s", broad_peaks_file))
broad_peaks <- readRDS(broad_peaks_file)
log_message(sprintf("Loading broad genes from %s", broad_genes_file))
broad_genes <- read.csv(broad_genes_file)

broad_viz <- create_visualizations(output_dir, "broad", broad_peaks, broad_genes)

log_message("Visualization process completed successfully")