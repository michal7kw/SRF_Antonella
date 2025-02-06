#!/usr/bin/env Rscript

#####################################################################
# Downstream Analysis Script for V5 and H2AK119Ub ChIP-seq Data
#####################################################################

# DESCRIPTION:
# This script performs downstream analysis of the cross-referenced ChIP-seq data,
# including genomic annotation, TSS distribution analysis, and gene list comparisons.
# It generates comprehensive visualizations and reports to characterize the binding
# patterns and genomic context of V5 and H2AK119Ub peaks.

# INPUT FILES:
# - Processed peak categories from cross-reference analysis:
#   * results/narrow_broad_categorized_peaks.rds
#   * results/broad_broad_categorized_peaks.rds
#
# Each RDS file contains GRanges objects with the following peak categories:
# - v5_with_h2a: V5 peaks overlapping H2AK119Ub
# - v5_only: V5 peaks without H2AK119Ub
# - h2a_with_v5: H2AK119Ub peaks overlapping V5
# - h2a_only: H2AK119Ub peaks without V5

# OUTPUT FILES:
# - Genomic annotation plots:
#   * results/downstream_analysis/*_genomic_annotation.pdf
#   * results/downstream_analysis/*_comparison.pdf
#
# - Gene lists:
#   * results/downstream_analysis/*_common_genes.txt
#   * results/downstream_analysis/*_v5_with_h2a_specific_genes.txt
#   * results/downstream_analysis/*_v5_only_specific_genes.txt
#
# - Analysis summary:
#   * results/downstream_analysis/analysis_summary.txt

# Load required libraries for genomic analysis and visualization
suppressPackageStartupMessages({
    library(GenomicRanges)    # For handling genomic intervals
    library(ChIPseeker)       # For ChIP peak annotation
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Gene annotations for hg38
    library(org.Hs.eg.db)     # Gene ID mappings for human
    library(ggplot2)          # For creating plots
    library(ComplexHeatmap)   # For heatmap visualization
    library(circlize)         # For circular plots
    library(GenomeInfoDb)     # For genome information
    library(rtracklayer)      # For importing genomic files
    library(RColorBrewer)     # For color palettes
})

# Define input/output directories
input_dir <- "results"
output_dir <- file.path(input_dir, "downstream_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load gene annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Set global theme for consistent plot appearance
theme_set(theme_minimal(base_size = 12) +
          theme(text = element_text(family = "sans")))

# Function to convert Entrez IDs to gene symbols with robust error handling
entrez_to_symbol <- function(entrez_ids) {
    # Remove any NA or empty values
    entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]
    
    if(length(entrez_ids) == 0) return(character(0))
    
    # Convert to symbols with error handling
    symbols <- tryCatch({
        mapIds(org.Hs.eg.db,
               keys = entrez_ids,
               column = "SYMBOL",
               keytype = "ENTREZID",
               multiVals = "first")
    }, error = function(e) {
        warning("Error converting Entrez IDs to symbols: ", e$message)
        return(NULL)
    })
    
    if(is.null(symbols)) return(entrez_ids)
    
    # Replace NAs with original IDs
    symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
    return(symbols)
}

# Custom plotting wrapper for ChIPseeker plots with enhanced visualization
custom_plot_annotation <- function(peak_annot, title, plot_type = "pie") {
    # Extract the count part
    count_text <- sprintf("(n = %d peaks)", length(peak_annot@anno))
    # Clean the title from any existing count
    title <- gsub("\\s*\\(n = \\d+ peaks\\)", "", title)
    
    # Set up PDF device with better resolution
    if(plot_type == "pie") {
        # Create color palette for consistent visualization
        colors <- c(
            "Promoter (<=1kb)" = "#2166AC",
            "Promoter (1-2kb)" = "#4393C3",
            "Promoter (2-3kb)" = "#92C5DE",
            "5' UTR" = "#D1E5F0",
            "3' UTR" = "#F7F7F7",
            "1st Exon" = "#FDDBC7",
            "Other Exon" = "#F4A582",
            "1st Intron" = "#D6604D",
            "Other Intron" = "#B2182B",
            "Downstream (<=300)" = "#67001F",
            "Distal Intergenic" = "#053061"
        )
        
        # Prepare data for plotting
        plot_data <- data.frame(
            feature = as.character(peak_annot@annoStat$Feature),
            frequency = peak_annot@annoStat$Frequency
        )
        
        # Create pie chart using ggplot2 with enhanced aesthetics
        p <- ggplot(plot_data, aes(x = "", y = frequency, fill = feature)) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar("y", start = 0) +
            scale_fill_manual(values = colors, name = "Feature") +
            labs(title = title,
                 subtitle = count_text,
                 x = NULL, y = NULL) +
            theme_minimal(base_size = 12) +
            theme(
                text = element_text(family = "sans"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
                legend.position = "right",
                legend.text = element_text(size = 9),
                panel.grid = element_blank(),
                axis.text = element_blank()
            )
        
        print(p)
        
    } else if(plot_type == "tss") {
        # Plot TSS distribution with improved visualization parameters
        p <- plotDistToTSS(peak_annot,
                        xlab = "Distance to TSS (bp)",
                        ylab = "Percentage of Peaks (%)",
                        title = title) +
            theme_minimal(base_size = 12) +
            theme(
                text = element_text(family = "sans"),
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 10),
                panel.grid.minor = element_line(color = "gray90", linetype = "dotted"),
                panel.grid.major = element_line(color = "gray85")
            )
        
        # Add subtitle with peak count
        p <- p + labs(subtitle = count_text)
        
        print(p)
    }
}

# Function to perform comprehensive peak region annotation analysis
analyze_peak_regions <- function(peaks, name, title_prefix) {
    if(length(peaks) == 0) {
        warning(sprintf("No peaks provided for analysis: %s", name))
        return(NULL)
    }
    
    # Annotate peaks using ChIPseeker
    peak_annot <- annotatePeak(peaks, 
                              TxDb = txdb,
                              annoDb = "org.Hs.eg.db",
                              level = "gene")
    
    # Save annotation results to PDF
    pdf(file.path(output_dir, sprintf("%s_genomic_annotation.pdf", name)),
        width = 10, height = 7)
    
    # Generate genomic annotation pie chart
    custom_plot_annotation(peak_annot, 
                         sprintf("%s", title_prefix),
                         "pie")
    
    # Generate TSS distribution plot
    custom_plot_annotation(peak_annot,
                         sprintf("%s\nDistance to TSS Distribution", title_prefix),
                         "tss")
    
    dev.off()
    
    return(peak_annot)
}

# Main analysis function coordinating all analysis steps
main <- function() {
    # Load processed peak categories from cross-reference analysis
    categories <- list(
        narrow_broad = readRDS(file.path(input_dir, "narrow_broad_categorized_peaks.rds")),
        broad_broad = readRDS(file.path(input_dir, "broad_broad_categorized_peaks.rds"))
    )
    
    # Analyze each peak category
    for (cat_name in names(categories)) {
        cat_peaks <- categories[[cat_name]]
        
        # Create readable titles for plots
        cat_title <- gsub("_", " ", toupper(cat_name))
        
        # Analyze V5 peaks with H2AK119Ub overlap
        v5_with_h2a_results <- analyze_peak_regions(
            cat_peaks$v5_with_h2a,
            sprintf("%s_v5_with_h2a", cat_name),
            sprintf("%s\nV5 Peaks with H2AK119Ub", cat_title)
        )
        
        # Analyze V5-only peaks
        v5_only_results <- analyze_peak_regions(
            cat_peaks$v5_only,
            sprintf("%s_v5_only", cat_name),
            sprintf("%s\nV5-only Peaks", cat_title)
        )
        
        # Compare gene lists to identify unique and shared targets
        common_genes <- intersect(v5_with_h2a_results@anno$geneId, v5_only_results@anno$geneId)
        v5_with_h2a_specific <- setdiff(v5_with_h2a_results@anno$geneId, v5_only_results@anno$geneId)
        v5_only_specific <- setdiff(v5_only_results@anno$geneId, v5_with_h2a_results@anno$geneId)
        
        # Convert Entrez IDs to more readable gene symbols
        common_genes_symbols <- entrez_to_symbol(common_genes)
        v5_with_h2a_specific_symbols <- entrez_to_symbol(v5_with_h2a_specific)
        v5_only_specific_symbols <- entrez_to_symbol(v5_only_specific)
        
        # Save gene lists to text files
        write.table(common_genes_symbols,
                   file = file.path(output_dir, sprintf("%s_common_genes.txt", cat_name)),
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(v5_with_h2a_specific_symbols,
                   file = file.path(output_dir, sprintf("%s_v5_with_h2a_specific_genes.txt", cat_name)),
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(v5_only_specific_symbols,
                   file = file.path(output_dir, sprintf("%s_v5_only_specific_genes.txt", cat_name)),
                   row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        # Generate comparison plots
        plot_comparison <- function(with_h2a_annot, only_annot, name) {
            pdf(file.path(output_dir, sprintf("%s_comparison.pdf", name)))
            
            # Compare genomic distributions using pie charts
            par(mfrow = c(2,1), oma = c(0,0,2,0))
            
            # Set up layout for two pie charts with legends
            layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths = c(6,4))
            
            # First pie chart with legend (V5 with H2AK119Ub)
            pie_data1 <- with_h2a_annot@annoStat$Frequency
            labels1 <- with_h2a_annot@annoStat$Feature
            colors1 <- rainbow(length(pie_data1))
            count_text1 <- sprintf("(n = %d peaks)", length(with_h2a_annot@anno))
            
            par(mar = c(4,2,4,1))
            pie(pie_data1, labels = NA, 
                main = sprintf("%s\nV5 Peaks with H2AK119Ub",
                              cat_title),
                col = colors1,
                cex.main = 0.9)
            mtext(count_text1, side = 1, line = 1, cex = 0.8)
            
            par(mar = c(4,0,4,2))
            plot.new()
            legend("center", 
                   legend = paste0(labels1, " (", round(pie_data1, 1), "%)"),
                   fill = colors1,
                   bty = "n")
            
            # Second pie chart with legend (V5-only)
            pie_data2 <- only_annot@annoStat$Frequency
            labels2 <- only_annot@annoStat$Feature
            colors2 <- rainbow(length(pie_data2))
            count_text2 <- sprintf("(n = %d peaks)", length(only_annot@anno))
            
            par(mar = c(4,2,4,1))
            pie(pie_data2, labels = NA,
                main = sprintf("%s\nV5-only Peaks",
                              cat_title),
                col = colors2,
                cex.main = 0.9)
            mtext(count_text2, side = 1, line = 1, cex = 0.8)
            
            par(mar = c(4,0,4,2))
            plot.new()
            legend("center", 
                   legend = paste0(labels2, " (", round(pie_data2, 1), "%)"),
                   fill = colors2,
                   bty = "n")
            
            # Reset layout for TSS distribution plots
            layout(1)
            par(mfrow = c(2,1), mar = c(5,4,6,2))
            
            # Compare TSS distributions
            plotDistToTSS(with_h2a_annot, 
                          title = sprintf("%s\nV5 Peaks with H2AK119Ub",
                                cat_title))
            mtext(count_text1, side = 3, line = 1, cex = 0.8)
            
            plotDistToTSS(only_annot, 
                          title = sprintf("%s\nV5-only Peaks",
                                cat_title))
            mtext(count_text2, side = 3, line = 1, cex = 0.8)
            
            dev.off()
        }
        
        # Generate comparison plots for each category
        plot_comparison(v5_with_h2a_results,
                       v5_only_results,
                       sprintf("%s_annotation", cat_name))
    }
    
    # Generate comprehensive analysis summary report
    sink(file.path(output_dir, "analysis_summary.txt"))
    cat("Downstream Analysis Summary\n")
    cat("=========================\n\n")
    
    for (cat_name in names(categories)) {
        cat_title <- gsub("_", " ", toupper(cat_name))
        cat(sprintf("\n%s Analysis:\n", cat_title))
        cat(sprintf("-----------------%s\n", paste(rep("-", nchar(cat_title)), collapse="")))
        
        cat_peaks <- categories[[cat_name]]
        
        cat(sprintf("V5 peaks with H2AK119Ub: %d\n", length(cat_peaks$v5_with_h2a)))
        cat(sprintf("V5-only peaks: %d\n", length(cat_peaks$v5_only)))
        cat(sprintf("H2AK119Ub peaks with V5: %d\n", length(cat_peaks$h2a_with_v5)))
        cat(sprintf("H2AK119Ub-only peaks: %d\n", length(cat_peaks$h2a_only)))
    }
    
    sink()
}

# Execute main analysis function
main() 