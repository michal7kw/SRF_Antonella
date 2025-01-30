source("scripts/utils.R")

#' Create comprehensive visualizations for peak analysis
#' @param peak_type Character, either "broad" or "narrow"
#' @param diff_peaks DiffBind results
#' @param gene_list Gene list from annotation
create_visualizations <- function(peak_type, diff_peaks, gene_list) {
    log_message(sprintf("Creating visualizations for %s peaks", peak_type))
    
    # Define base output directory
    base_dir <- file.path("analysis", paste0("plots_", peak_type))
    
    # Create output directories with full paths
    dirs <- c(
        base_dir,
        file.path(base_dir, "peak_analysis"),
        file.path(base_dir, "gene_analysis"),
        file.path(base_dir, "summary_statistics")
    )
    sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
    
    # Convert peaks to data frame if needed
    if (is(diff_peaks, "GRanges")) {
        diff_peaks <- as.data.frame(diff_peaks)
    }
    
    # 1. Peak Distribution Plot
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
             file.path(base_dir, "peak_analysis", "peak_distribution.pdf"))
    
    # 2. Enhanced MA Plot
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
             file.path(base_dir, "peak_analysis", "ma_plot_enhanced.pdf"))
    
    # 3. Enhanced Volcano Plot
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
             file.path(base_dir, "peak_analysis", "volcano_plot_enhanced.pdf"))
    
    # 4. Gene Distance to TSS Distribution
    log_message("Creating TSS distance distribution plot...")
    tss_dist_plot <- ggplot(gene_list, aes(x = distanceToTSS)) +
        geom_histogram(bins = 50, fill = "darkred", color = "black", alpha = 0.7) +
        theme_minimal() +
        labs(x = "Distance to TSS (bp)",
             y = "Frequency",
             title = paste("Distribution of Peak Distances to TSS -", peak_type)) +
        scale_x_continuous(labels = scales::comma)
    
    save_plot(tss_dist_plot,
             file.path(base_dir, "gene_analysis", "tss_distribution.pdf"))
    
    # 5. Create summary statistics
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
              file.path(base_dir, "summary_statistics", "summary_stats.csv"),
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

# Process both peak types
log_message("Starting visualization process...")

# Process broad peaks
broad_peaks <- readRDS("analysis/diffbind_broad/significant_peaks.rds")
broad_genes <- read.csv("analysis/gene_lists_broad/YAF_enriched_genes_broad_full.csv")
broad_viz <- create_visualizations("broad", broad_peaks, broad_genes)

# # Process narrow peaks
# narrow_peaks <- readRDS("analysis/diffbind_narrow/significant_peaks.rds")
# narrow_genes <- read.csv("analysis/gene_lists_narrow/YAF_enriched_genes_narrow_full.csv")
# narrow_viz <- create_visualizations("narrow", narrow_peaks, narrow_genes)

# # Create combined visualization report
# log_message("Creating combined visualization report...")
# pdf("analysis/plots_combined/combined_visualization_report.pdf", width = 12, height = 8)

# # Compare broad vs narrow peaks
# grid.arrange(
#     broad_viz$peak_dist_plot, narrow_viz$peak_dist_plot,
#     broad_viz$volcano_plot, narrow_viz$volcano_plot,
#     ncol = 2,
#     top = "Comparison of Broad and Narrow Peak Analysis"
# )

# dev.off()

log_message("Visualization process completed successfully") 