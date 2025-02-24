# This script performs annotation and enrichment analysis on peaks identified from
# differential binding analysis. It takes the peaks, annotates them relative to genes,
# performs GO enrichment analysis, and generates various visualizations.
#
# Input files:
# - analysis/7_differential_binding/diffbind_broad/all_peaks.rds: GRanges object with all peaks from DiffBind analysis
# - analysis/7_differential_binding/diffbind_broad/significant_peaks.rds: GRanges object with significant differential peaks
#
# Output files:
# In analysis/8_annotation_and_enrichment/annotation_broad/:
#   figures/
#     - annotation_plots.pdf: Peak annotation visualizations (pie chart, TSS distance)
#     - go_enrichment_plots.pdf: GO term enrichment plots (dotplot, emap, cnet)
#   tables/
#     - peak_annotation.csv: Detailed peak annotations
#     - go_enrichment.csv: GO enrichment analysis results
#   peak_annotation.rds: R object with full annotation data
#
# In analysis/8_annotation_and_enrichment/gene_lists_broad/:
#   - YAF_enriched_genes_broad_full.csv: All enriched genes with details
#   - YAF_enriched_genes_broad_symbols.txt: List of gene symbols only
#
# Dependencies:
# - ChIPseeker for peak annotation
# - clusterProfiler for GO enrichment
# - org.Hs.eg.db for gene ID mapping
# - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
# - ggupset (optional) for upset plots

# Load necessary libraries
suppressPackageStartupMessages({
    library(ChIPseeker)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(tidyverse)
    library(GenomicRanges)
})

source("scripts/utils.R")

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
    all_peaks_file <- file.path("analysis", "7_differential_binding", paste0("diffbind_", peak_type), "all_peaks.rds")
    significant_peaks_file <- file.path("analysis", "7_differential_binding", paste0("diffbind_", peak_type), "significant_peaks.rds")

    # Load peaks
    log_message(sprintf("Loading peaks from %s and %s", all_peaks_file, significant_peaks_file))
    all_peaks <- readRDS(all_peaks_file)
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
    all_peaks_gr <- standardize_chromosomes(all_peaks)
    significant_peaks_gr <- standardize_chromosomes(significant_peaks)

    # Annotate peaks relative to genomic features (promoters, UTRs, exons, etc)
    log_message("Performing peak annotation...")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno <- annotatePeak(significant_peaks_gr,
        TxDb = txdb,
        tssRegion = c(-3000, 3000), # Define promoter region
        verbose = FALSE
    )

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

    annotation_plots_pdf <- file.path(figures_dir, "annotation_plots.pdf")
    pdf(annotation_plots_pdf, width = 10, height = 8)
    print(plots$anno_pie)
    print(plots$dist_tss)
    if (!is.null(plots$upset)) {
        print(plots$upset)
    }
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
    significant_peaks_file <- file.path("analysis", "7_differential_binding", "diffbind_broad", "significant_peaks.rds")
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