# This script performs annotation and enrichment analysis on peaks identified from
# differential binding analysis. It takes the peaks, annotates them relative to genes,
# performs GO enrichment analysis, and generates various visualizations.
#
# Input files:
# - analysis/diffbind_broad/significant_peaks.rds: GRanges object with differential peaks
# - analysis/diffbind_narrow/significant_peaks.rds (optional): Narrow peaks if analyzing both types
#
# Output files:
# In analysis/annotation_{peak_type}/:
#   figures/
#     - annotation_plots.pdf: Peak annotation visualizations (pie chart, TSS distance)
#     - go_enrichment_plots.pdf: GO term enrichment plots (dotplot, emap, cnet)
#   tables/
#     - peak_annotation.csv: Detailed peak annotations
#     - go_enrichment.csv: GO enrichment analysis results
#   peak_annotation.rds: R object with full annotation data
#
# In analysis/gene_lists_{peak_type}/:
#   - YAF_enriched_genes_{peak_type}_full.csv: All enriched genes with details
#   - YAF_enriched_genes_{peak_type}_symbols.txt: List of gene symbols only
#
# Dependencies:
# - ChIPseeker for peak annotation
# - clusterProfiler for GO enrichment
# - org.Hs.eg.db for gene ID mapping
# - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
# - ggupset (optional) for upset plots

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
#' @param peak_type Character, either "broad" or "narrow"
#' @param peaks GRanges object from differential binding analysis
annotate_and_enrich <- function(peak_type = c("broad", "narrow"), peaks) {
    peak_type <- match.arg(peak_type)
    log_message(sprintf("Starting annotation and enrichment analysis for %s peaks", peak_type))
    
    # Create output directories
    dirs <- c(
        file.path("analysis", paste0("annotation_", peak_type), "figures"),
        file.path("analysis", paste0("annotation_", peak_type), "tables"),
        file.path("analysis", paste0("gene_lists_", peak_type))
    )
    create_dirs(dirs)
    
    # Standardize chromosome names to UCSC style (e.g., chr1 instead of 1)
    peaks_gr <- standardize_chromosomes(peaks)
    
    # Annotate peaks relative to genomic features (promoters, UTRs, exons, etc)
    log_message("Performing peak annotation...")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno <- annotatePeak(peaks_gr,
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),  # Define promoter region
                            verbose = FALSE)
    
    # Save annotation results
    saveRDS(peak_anno, 
            file.path("analysis", paste0("annotation_", peak_type), "peak_annotation.rds"))
    write.csv(as.data.frame(peak_anno),
             file.path("analysis", paste0("annotation_", peak_type), "tables", "peak_annotation.csv"),
             row.names = FALSE)
    
    # Generate visualization plots for peak annotations
    log_message("Generating annotation plots...")
    plots <- list(
        anno_pie = plotAnnoPie(peak_anno),  # Distribution of peak locations
        dist_tss = plotDistToTSS(peak_anno) # Distance of peaks to nearest TSS
    )
    
    # Create upset plot if package is available
    if (requireNamespace("ggupset", quietly = TRUE)) {
        plots$upset <- upsetplot(peak_anno, vennpie=TRUE)
    } else {
        log_message("Package 'ggupset' not available. Skipping upset plot.", level="WARNING")
    }
    
    pdf_file <- file.path("analysis", paste0("annotation_", peak_type), "figures", "annotation_plots.pdf")
    save_plot(plots$anno_pie, pdf_file, width=10, height=8)
    
    # Extract genes associated with YAF-enriched peaks
    log_message("Extracting YAF-enriched genes...")
    genes_df <- as.data.frame(peak_anno)
    
    # Convert Entrez IDs to gene symbols for better readability
    gene_ids <- genes_df$geneId[!is.na(genes_df$geneId)]
    gene_symbols <- mapIds(org.Hs.eg.db,
                          keys = gene_ids,
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first")
    
    # Create comprehensive gene list with annotations
    valid_indices <- which(!is.na(genes_df$geneId))
    gene_list <- data.frame(
        ENTREZID = gene_ids,
        SYMBOL = gene_symbols,
        distanceToTSS = genes_df$distanceToTSS[valid_indices],
        annotation = genes_df$annotation[valid_indices],
        fold_change = peaks_gr$Fold[valid_indices],
        stringsAsFactors = FALSE
    )
    
    # Clean up and sort gene list by fold change
    gene_list <- gene_list[!is.na(gene_list$SYMBOL), ]
    gene_list <- gene_list[order(-gene_list$fold_change), ]
    
    # Save gene lists in different formats
    write.csv(gene_list,
              file.path("analysis", paste0("gene_lists_", peak_type), 
                       paste0("YAF_enriched_genes_", peak_type, "_full.csv")),
              row.names = FALSE)
    
    writeLines(unique(gene_list$SYMBOL),
              file.path("analysis", paste0("gene_lists_", peak_type),
                       paste0("YAF_enriched_genes_", peak_type, "_symbols.txt")))
    
    # Perform GO enrichment analysis on YAF-enriched genes
    log_message("Performing GO enrichment analysis...")
    ego <- enrichGO(gene = unique(gene_list$ENTREZID),
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",  # Biological Process ontology
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
    
    if (!is.null(ego) && nrow(ego) > 0) {
        # Save GO enrichment results
        write.csv(as.data.frame(ego),
                 file.path("analysis", paste0("annotation_", peak_type), 
                          "tables", "go_enrichment.csv"),
                 row.names = FALSE)
        
        # Generate various GO enrichment visualizations
        pdf(file.path("analysis", paste0("annotation_", peak_type), 
                     "figures", "go_enrichment_plots.pdf"))
        
        # Dotplot showing enriched terms
        print(dotplot(ego, showCategory = 20))
        
        # Network plot showing term similarity
        tryCatch({
            ego <- pairwise_termsim(ego)
            print(emapplot(ego, showCategory = 50))
        }, error = function(e) {
            log_message("Could not create emapplot. Skipping...", level="WARNING")
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

# Process both peak types if available
log_message("Loading differential binding results...")

# Process broad peaks
broad_peaks <- readRDS("analysis/diffbind_broad/significant_peaks.rds")
broad_results <- annotate_and_enrich("broad", broad_peaks)

# Process narrow peaks (commented out as not currently used)
# narrow_peaks <- readRDS("analysis/diffbind_narrow/significant_peaks.rds")
# narrow_results <- annotate_and_enrich("narrow", narrow_peaks)

# Create summary comparing broad vs narrow peaks (commented out)
# log_message("Creating combined analysis summary...")
# combined_summary <- data.frame(
#     Analysis_Type = c("Broad", "Narrow"),
#     Total_Peaks = c(length(broad_peaks), length(narrow_peaks)),
#     Unique_Genes = c(length(unique(broad_results$gene_list$SYMBOL)),
#                     length(unique(narrow_results$gene_list$SYMBOL))),
#     GO_Terms = c(nrow(broad_results$go_enrichment), 
#                 nrow(narrow_results$go_enrichment))
# )

# write.csv(combined_summary,
#           "analysis/annotation_combined/combined_summary.csv",
#           row.names = FALSE)

log_message("Analysis completed successfully") 