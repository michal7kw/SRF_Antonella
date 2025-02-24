#!/usr/bin/env Rscript

#####################################################################
# Downstream Analysis Script for V5 and H2AK119Ub ChIP-seq Data
#####################################################################

# DESCRIPTION:
# This script performs downstream analysis of the cross-referenced ChIP-seq data,
# including genomic annotation, TSS distribution analysis, and gene list comparisons.
#
# INPUT FILES:
# - Processed peak categories from cross-reference analysis:
#   * results/broad_broad_categorized_peaks.rds
#
# Each RDS file contains GRanges objects with the following peak categories:
# - v5_with_h2a: V5 peaks overlapping H2AK119Ub
# - v5_only: V5 peaks without H2AK119Ub
# - h2a_with_v5: H2AK119Ub peaks overlapping V5
# - h2a_only: H2AK119Ub peaks without V5
#
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
# - CSV files:
#   * results/downstream_analysis/*_v5_with_h2a_overlapping_peaks.csv
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
        # Wrap TSS plotting in tryCatch to gracefully handle errors
        tryCatch({
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
        }, error = function(e) {
            warning("Error creating TSS distribution plot: ", e$message)
        })
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

# Main analysis function coordinating all processing steps
main_analysis <- function() {
    # Load processed peak categories from cross-reference analysis
    categorized_peaks <- readRDS(file.path(input_dir, "broad_broad_categorized_peaks.rds"))
    
    # Analyze each peak category
    categories <- list(
        v5_with_h2a = list(
            peaks = categorized_peaks$v5_with_h2a,
            name = "v5_with_h2a",
            title = "V5 Peaks with H2AK119Ub"
        ),
        v5_only = list(
            peaks = categorized_peaks$v5_only,
            name = "v5_only",
            title = "V5-only Peaks"
        ),
        h2a_with_v5 = list(
            peaks = categorized_peaks$h2a_with_v5,
            name = "h2a_with_v5",
            title = "H2AK119Ub Peaks with V5"
        ),
        h2a_only = list(
            peaks = categorized_peaks$h2a_only,
            name = "h2a_only",
            title = "H2AK119Ub-only Peaks"
        )
    )
    
    # Process each category
    results <- list()
    for (cat_name in names(categories)) {
        cat_info <- categories[[cat_name]]
        results[[cat_name]] <- analyze_peak_regions(
            cat_info$peaks,
            cat_info$name,
            cat_info$title
        )
    }
    
    # Extract gene lists
    for (cat_name in names(results)) {
        if (!is.null(results[[cat_name]])) {
            genes <- unique(results[[cat_name]]@anno$geneId)
            gene_symbols <- entrez_to_symbol(genes)
            
            # Save gene lists
            write.table(
                gene_symbols,
                file = file.path(output_dir, sprintf("%s_genes.txt", cat_name)),
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )
            
            # Save detailed peak information
            if (cat_name == "v5_with_h2a") {
                peak_info <- as.data.frame(results[[cat_name]]@anno)
                peak_info$geneSymbol <- entrez_to_symbol(peak_info$geneId)
                write.csv(
                    peak_info,
                    file = file.path(output_dir, "v5_with_h2a_overlapping_peaks.csv"),
                    row.names = FALSE
                )
            }
        }
    }
    
    # Generate comparison plots
    pdf(file.path(output_dir, "peak_comparisons.pdf"), width = 12, height = 8)
    
    # Create comparison plots between categories
    plot_category_comparison(results)
    
    dev.off()
    
    # Generate analysis summary
    generate_analysis_summary(results, categorized_peaks)
}

# Function to plot category comparisons
plot_category_comparison <- function(results) {
    # Extract peak counts
    peak_counts <- sapply(results, function(x) {
        if (is.null(x)) return(0)
        length(x@anno)
    })
    
    # Create barplot of peak counts
    barplot_data <- data.frame(
        Category = names(peak_counts),
        Count = unname(peak_counts)
    )
    
    p1 <- ggplot(barplot_data, aes(x = Category, y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Peak Counts by Category",
             x = "Category",
             y = "Number of Peaks")
    
    print(p1)
    
    # Create genomic distribution comparison
    if (!is.null(results$v5_with_h2a) && !is.null(results$v5_only)) {
        genomic_dist_data <- rbind(
            data.frame(
                Category = "V5 with H2AK119Ub",
                Feature = results$v5_with_h2a@annoStat$Feature,
                Frequency = results$v5_with_h2a@annoStat$Frequency
            ),
            data.frame(
                Category = "V5 only",
                Feature = results$v5_only@annoStat$Feature,
                Frequency = results$v5_only@annoStat$Frequency
            )
        )
        
        p2 <- ggplot(genomic_dist_data, aes(x = Feature, y = Frequency, fill = Category)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "Genomic Distribution Comparison",
                 x = "Genomic Feature",
                 y = "Frequency (%)")
        
        print(p2)
    }
}

# Function to generate analysis summary
generate_analysis_summary <- function(results, categorized_peaks) {
    sink(file.path(output_dir, "analysis_summary.txt"))
    
    cat("Downstream Analysis Summary\n")
    cat("==========================\n\n")
    
    cat("Peak Counts:\n")
    cat("-----------\n")
    for (cat_name in names(categorized_peaks)) {
        cat(sprintf("%s: %d peaks\n", 
                   gsub("_", " ", toupper(cat_name)), 
                   length(categorized_peaks[[cat_name]])))
    }
    
    cat("\nGenomic Distribution Summary:\n")
    cat("--------------------------\n")
    for (cat_name in names(results)) {
        if (!is.null(results[[cat_name]])) {
            cat(sprintf("\n%s:\n", gsub("_", " ", toupper(cat_name))))
            print(results[[cat_name]]@annoStat)
        }
    }
    
    cat("\nAnalysis Details:\n")
    cat("----------------\n")
    cat(sprintf("Analysis performed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Genome build: hg38\n"))
    cat(sprintf("Gene annotation database: %s\n", packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene")))
    
    sink()
}

# Execute main analysis
main_analysis() 