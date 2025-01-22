# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(DiffBind)
})

# Function to create plots for a specific peak type
create_plots <- function(peak_type) {
  # Create output directories
  plots_dir <- file.path("analysis", paste0("plots_", peak_type))
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read the DiffBind results with error handling
  diff_peaks_file <- file.path("analysis", paste0("diffbind_", peak_type), "differential_peaks.csv")
  
  tryCatch({
    # Read the CSV file
    diff_peaks <- read.csv(diff_peaks_file)
    
    # Add debugging information
    message(paste0("Processing ", peak_type, " peaks"))
    message("Total peaks loaded: ", nrow(diff_peaks))
    
    # Create MA plot with enhanced styling
    ma_plot <- ggplot(diff_peaks, aes(x = log2(Conc), y = Fold)) +
      geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
      scale_color_manual(values = c("grey50", "red3")) +
      scale_size_continuous(range = c(0.5, 2)) +
      labs(x = "log2 Mean Concentration",
           y = "log2 Fold Change (YAF/GFP)",
           title = paste("MA Plot: YAF vs GFP H2AK119Ub", peak_type, "Peaks"),
           subtitle = paste0("Significant peaks: ", sum(diff_peaks$FDR < 0.05))) +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    
    # Save MA plot
    ma_plot_file <- file.path(plots_dir, paste0("MA_plot_", peak_type, ".pdf"))
    ggsave(ma_plot_file, ma_plot, width = 8, height = 6)
    message("Created MA plot: ", ma_plot_file)
    
    # Create volcano plot
    volcano_plot <- ggplot(diff_peaks, aes(x = Fold, y = -log10(FDR))) +
      geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
      scale_color_manual(values = c("grey50", "red3")) +
      scale_size_continuous(range = c(0.5, 2)) +
      labs(x = "log2 Fold Change (YAF/GFP)",
           y = "-log10(FDR)",
           title = paste("Volcano Plot: YAF vs GFP H2AK119Ub", peak_type, "Peaks"),
           subtitle = paste0("Significant peaks: ", sum(diff_peaks$FDR < 0.05))) +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    
    # Save volcano plot
    volcano_plot_file <- file.path(plots_dir, paste0("volcano_plot_", peak_type, ".pdf"))
    ggsave(volcano_plot_file, volcano_plot, width = 8, height = 6)
    message("Created volcano plot: ", volcano_plot_file)
    
    # Calculate summary statistics
    summary_stats <- data.frame(
      Metric = c("Total Peaks", "Significant Peaks", "YAF-enriched", "GFP-enriched"),
      Count = c(
        nrow(diff_peaks),
        sum(diff_peaks$FDR < 0.05),
        sum(diff_peaks$FDR < 0.05 & diff_peaks$Fold > 0),
        sum(diff_peaks$FDR < 0.05 & diff_peaks$Fold < 0)
      )
    )
    
    # Save summary statistics
    summary_file <- file.path(plots_dir, paste0("summary_stats_", peak_type, ".txt"))
    write.table(summary_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message("Created summary statistics: ", summary_file)
    
    return(summary_stats)
    
  }, error = function(e) {
    stop(paste("Error processing", peak_type, "peaks:", e$message))
  })
}

# Process both narrow and broad peaks
narrow_stats <- create_plots("narrow")
broad_stats <- create_plots("broad")

# Create combined summary
combined_summary <- rbind(
  data.frame(Peak_Type = "Narrow", narrow_stats),
  data.frame(Peak_Type = "Broad", broad_stats)
)

# Save combined summary
combined_dir <- "analysis/plots_combined"
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
write.table(combined_summary, 
            file.path(combined_dir, "combined_summary_stats.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Analysis completed for both narrow and broad peaks")
