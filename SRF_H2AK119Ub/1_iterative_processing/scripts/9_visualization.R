# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(GenomicRanges)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)
})

# Create visualization directory if it doesn't exist
dir.create("results/visualization", recursive = TRUE, showWarnings = FALSE)

# Read differential binding results
diff_peaks <- readRDS("analysis/diffbind/significant_peaks.rds")
all_peaks <- readRDS("analysis/diffbind/all_peaks.rds")

# Create volcano plot with enhanced styling
pdf("results/visualization/volcano_plot.pdf", width = 8, height = 7)
ggplot(as.data.frame(all_peaks), aes(x = Fold, y = -log10(FDR))) +
    geom_point(aes(color = FDR < 0.05 & abs(Fold) > 1), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("gray70", "red3")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change (YAF/GFP)",
         y = "-log10(FDR)",
         title = "H2AK119ub Differential Binding",
         color = "Significant") +
    theme(legend.position = "bottom")
dev.off()

# Create peak distribution plot
peak_stats <- data.frame(
    Category = c("Total Peaks", "YAF-enriched", "GFP-enriched"),
    Count = c(
        nrow(all_peaks),
        sum(diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05),
        sum(diff_peaks$Fold < 0 & diff_peaks$FDR < 0.05)
    )
)

pdf("results/visualization/peak_distribution.pdf", width = 8, height = 6)
ggplot(peak_stats, aes(x = reorder(Category, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    theme_minimal() +
    labs(x = "Category",
         y = "Number of Peaks",
         title = "Distribution of H2AK119ub Peaks") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Read gene annotation results
yaf_genes <- read.table("analysis/gene_lists/YAF_enriched_H2AK119Ub_genes.txt",
                       header = TRUE, sep = "\t")

# Create distance to TSS distribution plot
pdf("results/visualization/tss_distribution.pdf", width = 8, height = 6)
ggplot(yaf_genes, aes(x = distance_to_TSS)) +
    geom_histogram(bins = 50, fill = "darkred", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(x = "Distance to TSS (bp)",
         y = "Frequency",
         title = "Distribution of Peak Distances to TSS") +
    scale_x_continuous(labels = scales::comma)
dev.off()

# Create summary statistics
summary_stats <- data.frame(
    Metric = c(
        "Total Peaks Analyzed",
        "Significant Peaks (FDR < 0.05)",
        "YAF-enriched Peaks",
        "GFP-enriched Peaks",
        "Unique Genes",
        "Median Distance to TSS"
    ),
    Value = c(
        nrow(all_peaks),
        nrow(diff_peaks),
        sum(diff_peaks$Fold > 0 & diff_peaks$FDR < 0.05),
        sum(diff_peaks$Fold < 0 & diff_peaks$FDR < 0.05),
        length(unique(yaf_genes$gene_symbol)),
        median(abs(yaf_genes$distance_to_TSS))
    )
)

write.csv(summary_stats,
          "results/visualization/analysis_summary.csv",
          row.names = FALSE)

# Read GO enrichment results if they exist
go_file <- "analysis/gene_lists/YAF_enriched_GO_terms.txt"
if (file.exists(go_file)) {
    go_terms <- read.table(go_file, header = TRUE, sep = "\t")
    
    # Create GO terms visualization
    top_go <- head(go_terms[order(go_terms$p.adjust), ], 15)
    
    pdf("results/visualization/top_GO_terms.pdf", width = 10, height = 8)
    ggplot(top_go, aes(x = reorder(Description, -log10(p.adjust)), 
                       y = -log10(p.adjust))) +
        geom_bar(stat = "identity", fill = "forestgreen", alpha = 0.7) +
        coord_flip() +
        theme_minimal() +
        labs(x = "GO Term",
             y = "-log10(Adjusted P-value)",
             title = "Top 15 Enriched GO Terms") +
        theme(axis.text.y = element_text(size = 8))
    dev.off()
}

print("Visualization completed successfully")