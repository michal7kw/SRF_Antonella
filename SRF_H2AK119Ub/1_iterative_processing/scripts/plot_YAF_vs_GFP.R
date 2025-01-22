# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(DiffBind)
})

# Create output directories
dir.create("analysis/plots", recursive = TRUE, showWarnings = FALSE)

# Read the DiffBind results with error handling
tryCatch({
  # Read the tab-delimited file instead of RDS
  diff_peaks <- read.table("analysis/diffbind/all_differential_peaks.txt", 
                          header = TRUE, sep = "\t")
  
  # Standardize chromosome names
  diff_peaks$Chr <- ifelse(!grepl("^chr", diff_peaks$Chr, ignore.case = TRUE),
                          paste0("chr", diff_peaks$Chr),
                          diff_peaks$Chr)
  diff_peaks$Chr <- sub("chrMT", "chrM", diff_peaks$Chr)
  diff_peaks$Chr <- sub("chr23", "chrX", diff_peaks$Chr)
  diff_peaks$Chr <- sub("chr24", "chrY", diff_peaks$Chr)
  
  message("Unique chromosome names: ", paste(unique(diff_peaks$Chr), collapse = ", "))
  message("Total peaks loaded: ", nrow(diff_peaks))
  
}, error = function(e) {
  stop("Error reading DiffBind results: ", e$message)
})

# Create MA plot with enhanced styling
ma_plot <- ggplot(diff_peaks, aes(x = log2(Conc), y = Fold)) +
  geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(x = "log2 Mean Concentration",
       y = "log2 Fold Change (YAF/GFP)",
       title = "MA Plot: YAF vs GFP H2AK119Ub Peaks",
       subtitle = paste0("Significant peaks: ", sum(diff_peaks$FDR < 0.05))) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

# Save MA plot with error handling
tryCatch({
  ggsave("analysis/plots/MA_plot_YAF_vs_GFP.pdf", ma_plot, width = 8, height = 6)
  message("MA plot saved successfully")
}, error = function(e) {
  warning("Error saving MA plot: ", e$message)
})

# Create enhanced volcano plot with better filtering
volcano_plot <- ggplot(diff_peaks, 
                      aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05 & abs(Fold) > 1,
                 size = abs(Fold)), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(x = "log2 Fold Change (YAF/GFP)",
       y = "-log10(FDR)",
       title = "Volcano Plot: YAF vs GFP H2AK119Ub Peaks",
       subtitle = paste0("Up: ", sum(diff_peaks$Fold > 1 & diff_peaks$FDR < 0.05),
                        " Down: ", sum(diff_peaks$Fold < -1 & diff_peaks$FDR < 0.05))) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

# Save volcano plot with error handling
tryCatch({
  ggsave("analysis/plots/volcano_plot_YAF_vs_GFP.pdf", volcano_plot, width = 8, height = 6)
  message("Volcano plot saved successfully")
}, error = function(e) {
  warning("Error saving volcano plot: ", e$message)
})

# Create profile plot if deepTools output exists
profile_file <- "analysis/visualization_merged/YAF_vs_GFP_matrix.gz"
if(file.exists(profile_file)) {
  tryCatch({
    profile_data <- read.table(profile_file, header = TRUE)
    
    if(nrow(profile_data) == 0) {
      warning("Profile data is empty")
    } else {
      profile_plot <- ggplot(profile_data, aes(x = Position, y = Signal, color = Condition)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = Signal - SEM, ymax = Signal + SEM, fill = Condition), 
                    alpha = 0.2) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        labs(x = "Distance from peak center (bp)",
             y = "Average H2AK119Ub signal",
             title = "H2AK119Ub Signal Profile",
             subtitle = "Averaged signal around peak centers") +
        theme_bw() +
        theme(
          legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10)
        )
      
      ggsave("analysis/plots/profile_plot_YAF_vs_GFP.pdf", profile_plot, width = 8, height = 6)
      message("Profile plot saved successfully")
    }
  }, error = function(e) {
    warning("Error creating profile plot: ", e$message)
  })
} else {
  message("Profile data file not found: ", profile_file)
}

# Save summary statistics
summary_stats <- data.frame(
  Metric = c("Total Peaks", "Significant Peaks", "YAF-enriched", "GFP-enriched"),
  Count = c(
    nrow(diff_peaks),
    sum(diff_peaks$FDR < 0.05),
    sum(diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05),
    sum(diff_peaks$Fold < 0 & diff_peaks$FDR < 0.05)
  )
)

write.table(summary_stats,
            "analysis/plots/peak_statistics.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

message("Summary statistics saved successfully")
print("Visualization completed successfully")
