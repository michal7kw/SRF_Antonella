# This script performs annotation and enrichment analysis on peaks identified from
# differential binding analysis. It takes the peaks, annotates them relative to genes,
# performs GO enrichment analysis, and generates various visualizations.
#
# Input files:
# - analysis/7_differential_binding/significant_peaks.rds: GRanges object with significant differential peaks
#
# Output files:
# In analysis/8_annotation_and_enrichment/annotation_broad/:
#   figures/
#     - annotation_plots.pdf: Peak annotation visualizations (pie chart, TSS distance)
#     - detailed_pie_chart.pdf: Detailed pie chart with genomic feature distribution
#     - go_enrichment_plots.pdf: GO term enrichment plots (dotplot, emap, cnet)
#   tables/
#     - peak_annotation.csv: Detailed peak annotations
#     - go_enrichment.csv: GO enrichment analysis results
#   peak_annotation.rds: R object with full annotation data
#
# In analysis/8_annotation_and_enrichment/gene_lists_broad/:
#   - YAF_enriched_genes_broad_full.csv: All enriched genes with details
#   - YAF_enriched_genes_broad_symbols.txt: List of gene symbols only
#   - YAF_enriched_genes_broad_promoters.txt: Genes associated with promoter regions
#   - YAF_enriched_genes_broad_promoters.csv: Genes associated with promoter regions (with details)
#   - YAF_enriched_genes_broad_promoters_1st_exon_intron.txt: Genes in promoters + 1st exon/intron
#   - YAF_enriched_genes_broad_distal_intergenic.txt: Genes associated with distal intergenic regions
#   - YAF_enriched_genes_broad_other_regions.txt: Genes in other genomic regions
#
# Dependencies:
# - ChIPseeker for peak annotation
# - clusterProfiler for GO enrichment
# - org.Hs.eg.db for gene ID mapping
# - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
# - ggupset (optional) for upset plots
# - ggplot2 for custom visualizations

# Load necessary libraries
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(tidyverse)
    library(GenomicRanges)
    library(ggplot2)
})

source("scripts/utils.R")

#' Create a detailed pie chart showing genomic feature distribution
#' @param anno_data Data frame with annotation data and fold change information
#' @param direction Optional parameter to filter by direction ("up", "down", or NULL for all)
#' @param title Optional title for the pie chart
#' @return ggplot object with detailed pie chart
create_detailed_pie_chart <- function(anno_data, direction = NULL, title = NULL) {
    # Make a copy of the annotation data
    anno_df <- anno_data
    
    # Filter by direction if specified
    if (!is.null(direction)) {
        if (direction == "up") {
            anno_df <- anno_df[anno_df$fold_change > 0, ]
        } else if (direction == "down") {
            anno_df <- anno_df[anno_df$fold_change < 0, ]
        }
    }
    
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
    anno_df$detailed_annotation[grepl("Exon", anno_df$annotation)] <- "Exon"
    anno_df$detailed_annotation[grepl("Intron", anno_df$annotation)] <- "Intron"
    
    # Downstream and intergenic
    anno_df$detailed_annotation[grepl("Downstream", anno_df$annotation)] <- "Downstream (<=300bp)"
    anno_df$detailed_annotation[grepl("Distal Intergenic", anno_df$annotation)] <- "Distal Intergenic"
    
    # Calculate percentages
    anno_summary <- anno_df %>%
        group_by(detailed_annotation) %>%
        summarise(count = n()) %>%
        mutate(percentage = count / sum(count) * 100) %>%
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
        "Downstream (<=300bp)" = "#6A3D9A",
        "Distal Intergenic" = "#FFFF99"
    )
    
    # Create the pie chart
    pie_chart <- ggplot(anno_summary, aes(x = "", y = percentage, fill = detailed_annotation)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = color_palette) +
        theme_void() +
        theme(legend.position = "right") +
        labs(fill = "Genomic Feature") +
        geom_text(aes(label = paste0(round(percentage, 2), "%")), 
                  position = position_stack(vjust = 0.5),
                  size = 3)
    
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
    # Define standard chromosomes (1-22, X, Y)
    standard_chroms <- paste0("chr", c(1:22, "X", "Y"))
    
    # Add chr prefix if missing
    seqlevels(gr) <- paste0("chr", sub("^chr", "", seqlevels(gr)))
    
    # Keep only standard chromosomes
    gr_filtered <- gr[seqnames(gr) %in% standard_chroms]
    
    # Update sequence levels
    seqlevels(gr_filtered) <- standard_chroms[standard_chroms %in% seqlevels(gr_filtered)]
    
    return(gr_filtered)
}

#' Perform peak annotation and enrichment analysis
annotate_and_enrich <- function(output_dir, peak_type = "broad") {
    log_message(sprintf("Starting annotation and enrichment analysis for %s peaks", peak_type))

    # Define input file paths
    significant_peaks_file <- file.path("analysis", "7_differential_binding", "significant_peaks.rds")

    # Load peaks
    significant_peaks <- readRDS(significant_peaks_file)

    # Define output directories
    annotation_dir <- file.path(output_dir, paste0("annotation_", peak_type))
    gene_lists_dir <- file.path(output_dir, paste0("gene_lists_", peak_type))
    figures_dir <- file.path(annotation_dir, "figures")
    tables_dir <- file.path(annotation_dir, "tables")

    # Create output directories
    dirs <- c(figures_dir, tables_dir, gene_lists_dir)
    create_dirs(dirs)

    # Standardize chromosome names to UCSC style (e.g., chr1 instead of 1)
    significant_peaks_gr <- standardize_chromosomes(significant_peaks)

    # Add fold change information to the annotation object for later filtering
    fold_change_values <- significant_peaks_gr$Fold

    # Annotate peaks relative to genomic features (promoters, UTRs, exons, etc)
    log_message("Performing peak annotation...")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno <- annotatePeak(significant_peaks_gr,
        TxDb = txdb,
        tssRegion = c(-3000, 3000), # Define promoter region
        verbose = FALSE
    )
    
    # Instead of modifying the peak_anno object directly, we'll store the fold change
    # in a separate data frame that we'll use for filtering in the pie chart function
    anno_df <- as.data.frame(peak_anno)
    anno_df$fold_change <- fold_change_values
    
    # Save the modified annotation data frame for later use
    fold_change_data <- anno_df

    # Save annotation results
    peak_annotation_file <- file.path(annotation_dir, "peak_annotation.rds")
    saveRDS(peak_anno, peak_annotation_file)

    peak_annotation_csv_file <- file.path(tables_dir, "peak_annotation.csv")
    write.csv(as.data.frame(peak_anno),
        peak_annotation_csv_file,
        row.names = FALSE
    )

    # Generate visualization plots for peak annotations
    log_message("Generating annotation plots...")
    plots <- list(
        anno_pie = plotAnnoPie(peak_anno), # Distribution of peak locations
        dist_tss = plotDistToTSS(peak_anno) # Distance of peaks to nearest TSS
    )

    # Create upset plot if package is available
    if (requireNamespace("ggupset", quietly = TRUE)) {
        plots$upset <- upsetplot(peak_anno, vennpie = TRUE)
    } else {
        log_message("Package 'ggupset' not available. Skipping upset plot.", level = "WARNING")
    }

    # Create detailed pie chart
    plots$detailed_pie <- create_detailed_pie_chart(fold_change_data, title = "All Peaks")
    
    # Create separate pie charts for up and down regulated genes
    plots$up_regulated_pie <- create_detailed_pie_chart(fold_change_data, direction = "up", 
                                                       title = "Up-regulated Peaks")
    plots$down_regulated_pie <- create_detailed_pie_chart(fold_change_data, direction = "down", 
                                                         title = "Down-regulated Peaks")

    annotation_plots_pdf <- file.path(figures_dir, "annotation_plots.pdf")
    pdf(annotation_plots_pdf, width = 10, height = 8)
    print(plots$anno_pie)
    print(plots$dist_tss)
    if (!is.null(plots$upset)) {
        print(plots$upset)
    }
    print(plots$detailed_pie)
    print(plots$up_regulated_pie)
    print(plots$down_regulated_pie)
    dev.off()

    # Save detailed pie charts to separate files for better visibility
    detailed_pie_pdf <- file.path(figures_dir, "detailed_pie_chart.pdf")
    pdf(detailed_pie_pdf, width = 10, height = 8)
    print(plots$detailed_pie)
    dev.off()
    
    # Save up-regulated pie chart
    up_regulated_pie_pdf <- file.path(figures_dir, "up_regulated_pie_chart.pdf")
    pdf(up_regulated_pie_pdf, width = 10, height = 8)
    print(plots$up_regulated_pie)
    dev.off()
    
    # Save down-regulated pie chart
    down_regulated_pie_pdf <- file.path(figures_dir, "down_regulated_pie_chart.pdf")
    pdf(down_regulated_pie_pdf, width = 10, height = 8)
    print(plots$down_regulated_pie)
    dev.off()

    # Extract genes associated with YAF-enriched peaks
    log_message("Extracting YAF-enriched genes...")
    genes_df <- as.data.frame(peak_anno)

    # Convert Entrez IDs to gene symbols for better readability
    gene_ids <- genes_df$geneId[!is.na(genes_df$geneId)]
    gene_symbols <- mapIds(org.Hs.eg.db,
        keys = gene_ids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )

    # Create comprehensive gene list with annotations
    valid_indices <- which(!is.na(genes_df$geneId))
    gene_list <- data.frame(
        ENTREZID = gene_ids,
        SYMBOL = gene_symbols,
        distanceToTSS = genes_df$distanceToTSS[valid_indices],
        annotation = genes_df$annotation[valid_indices],
        fold_change = significant_peaks_gr$Fold[valid_indices],
        p_value = significant_peaks_gr$`p-value`[valid_indices],
        FDR = significant_peaks_gr$FDR[valid_indices],
        stringsAsFactors = FALSE
    )

    # Clean up and sort gene list by fold change
    gene_list <- gene_list[!is.na(gene_list$SYMBOL), ]
    gene_list <- gene_list[order(-gene_list$fold_change), ]

    # Save gene lists in different formats
    yaf_enriched_genes_full_csv <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_full.csv"))
    write.csv(gene_list,
        yaf_enriched_genes_full_csv,
        row.names = FALSE
    )

    yaf_enriched_genes_symbols_txt <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_symbols.txt"))
    writeLines(unique(gene_list$SYMBOL),
        yaf_enriched_genes_symbols_txt
    )

    # Create categorized gene lists as requested
    log_message("Creating categorized gene lists...")
    
    # 1. Promoter regions
    promoter_genes <- gene_list[grepl("Promoter", gene_list$annotation), ]
    promoter_genes_file <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_promoters.txt"))
    writeLines(unique(promoter_genes$SYMBOL), promoter_genes_file)
    promoter_genes_csv <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_promoters.csv"))
    write.csv(promoter_genes, promoter_genes_csv, row.names = FALSE)
    
    # 2. Promoters + Exon/Intron
    promoter_exon_intron_genes <- gene_list[grepl("Promoter|Exon|Intron", gene_list$annotation), ]
    promoter_exon_intron_file <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_promoters_exons_introns.txt"))
    writeLines(unique(promoter_exon_intron_genes$SYMBOL), promoter_exon_intron_file)
    promoter_exon_intron_csv <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_promoters_exons_introns.csv"))
    write.csv(promoter_exon_intron_genes, promoter_exon_intron_csv, row.names = FALSE)
    
    # 3. Distal intergenic
    distal_intergenic_genes <- gene_list[grepl("Distal Intergenic", gene_list$annotation), ]
    distal_intergenic_file <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_distal_intergenic.txt"))
    writeLines(unique(distal_intergenic_genes$SYMBOL), distal_intergenic_file)
    distal_intergenic_csv <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_distal_intergenic.csv"))
    write.csv(distal_intergenic_genes, distal_intergenic_csv, row.names = FALSE)
    
    # 4. Other regions (not in the above categories)
    other_genes <- gene_list[!grepl("Promoter|Exon|Intron|Distal Intergenic", gene_list$annotation), ]
    other_genes_file <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_other_regions.txt"))
    writeLines(unique(other_genes$SYMBOL), other_genes_file)
    other_genes_csv <- file.path(gene_lists_dir, paste0("YAF_enriched_genes_", peak_type, "_other_regions.csv"))
    write.csv(other_genes, other_genes_csv, row.names = FALSE)

    # Perform GO enrichment analysis on YAF-enriched genes
    log_message("Performing GO enrichment analysis...")
    ego <- enrichGO(
        gene = unique(gene_list$ENTREZID),
        OrgDb = org.Hs.eg.db,
        ont = "BP", # Biological Process ontology
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )

    if (!is.null(ego) && nrow(ego) > 0) {
        # Save GO enrichment results
        go_enrichment_csv <- file.path(tables_dir, "go_enrichment.csv")
        write.csv(as.data.frame(ego),
            go_enrichment_csv,
            row.names = FALSE
        )

        # Generate various GO enrichment visualizations
        go_enrichment_plots_pdf <- file.path(figures_dir, "go_enrichment_plots.pdf")
        pdf(go_enrichment_plots_pdf)

        # Dotplot showing enriched terms
        print(dotplot(ego, showCategory = 20))

        # Network plot showing term similarity
        tryCatch({
            ego <- pairwise_termsim(ego)
            print(emapplot(ego, showCategory = 50))
        }, error = function(e) {
            log_message("Could not create emapplot. Skipping...", level = "WARNING")
        })

        # Gene-concept network plot
        print(cnetplot(ego, showCategory = 10))

        dev.off()
    }

    log_message(sprintf("Completed annotation and enrichment analysis for %s peaks", peak_type))
    return(list(
        annotation = peak_anno,
        gene_list = gene_list,
        go_enrichment = ego
    ))
}

# Main execution block
log_message("Starting annotation and enrichment pipeline...")

# Get the output directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Output directory must be provided as a command line argument.")
}
output_dir <- args[1]

# Perform annotation and enrichment analysis for broad peaks
broad_results <- annotate_and_enrich(output_dir, "broad")

# Create summary for broad peaks
log_message("Creating analysis summary...")
if (!is.null(broad_results)) {
    # Load the significant peaks to get the total number of peaks
    significant_peaks_file <- file.path("analysis", "7_differential_binding", "significant_peaks.rds")
    significant_peaks <- readRDS(significant_peaks_file)

    combined_summary <- data.frame(
        Analysis_Type = "Broad",
        Total_Peaks = length(significant_peaks),
        Unique_Genes = length(unique(broad_results$gene_list$SYMBOL)),
        GO_Terms = ifelse(!is.null(broad_results$go_enrichment), nrow(broad_results$go_enrichment), 0)
    )

    combined_summary_file <- file.path(output_dir, "combined_summary.csv")
    write.csv(combined_summary,
        combined_summary_file,
        row.names = FALSE
    )
} else {
    log_message("Skipping combined summary due to missing broad results.", level = "WARNING")
}

log_message("Analysis completed successfully")