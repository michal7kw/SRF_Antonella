source("scripts/utils.R")

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
    
    # Standardize chromosome names
    peaks_gr <- standardize_chromosomes(peaks)
    
    # Annotate peaks
    log_message("Performing peak annotation...")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno <- annotatePeak(peaks_gr,
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            verbose = FALSE)
    
    # Save annotation results
    saveRDS(peak_anno, 
            file.path("analysis", paste0("annotation_", peak_type), "peak_annotation.rds"))
    write.csv(as.data.frame(peak_anno),
             file.path("analysis", paste0("annotation_", peak_type), "tables", "peak_annotation.csv"),
             row.names = FALSE)
    
    # Create annotation plots
    log_message("Generating annotation plots...")
    plots <- list(
        anno_pie = plotAnnoPie(peak_anno),
        dist_tss = plotDistToTSS(peak_anno),
        upset = upsetplot(peak_anno, vennpie=TRUE)
    )
    
    pdf_file <- file.path("analysis", paste0("annotation_", peak_type), "figures", "annotation_plots.pdf")
    save_plot(plots$anno_pie, pdf_file, width=10, height=8)
    
    # Extract YAF-enriched genes
    log_message("Extracting YAF-enriched genes...")
    genes_df <- as.data.frame(peak_anno)
    
    # Map ENTREZID to SYMBOL
    gene_symbols <- mapIds(org.Hs.eg.db,
                          keys = genes_df$geneId,
                          column = "SYMBOL",
                          keytype = "ENTREZID",
                          multiVals = "first")
    
    # Create gene list
    gene_list <- data.frame(
        ENTREZID = genes_df$geneId,
        SYMBOL = gene_symbols,
        distanceToTSS = genes_df$distanceToTSS,
        annotation = genes_df$annotation,
        fold_change = peaks_gr$Fold,
        stringsAsFactors = FALSE
    ) %>%
        filter(!is.na(SYMBOL)) %>%
        arrange(desc(fold_change))
    
    # Save gene lists
    write.csv(gene_list,
              file.path("analysis", paste0("gene_lists_", peak_type), 
                       paste0("YAF_enriched_genes_", peak_type, "_full.csv")),
              row.names = FALSE)
    
    writeLines(unique(gene_list$SYMBOL),
              file.path("analysis", paste0("gene_lists_", peak_type),
                       paste0("YAF_enriched_genes_", peak_type, "_symbols.txt")))
    
    # Perform GO enrichment analysis
    log_message("Performing GO enrichment analysis...")
    ego <- enrichGO(gene = unique(gene_list$ENTREZID),
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
    
    if (!is.null(ego) && nrow(ego) > 0) {
        # Save GO results
        write.csv(as.data.frame(ego),
                 file.path("analysis", paste0("annotation_", peak_type), 
                          "tables", "go_enrichment.csv"),
                 row.names = FALSE)
        
        # Create GO plots
        pdf(file.path("analysis", paste0("annotation_", peak_type), 
                     "figures", "go_enrichment_plots.pdf"))
        print(dotplot(ego, showCategory = 20))
        print(emapplot(ego, showCategory = 50))
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

# Process both peak types
log_message("Loading differential binding results...")

# Process broad peaks
broad_peaks <- readRDS("analysis/diffbind_broad/significant_peaks.rds")
broad_results <- annotate_and_enrich("broad", broad_peaks)

# Process narrow peaks
narrow_peaks <- readRDS("analysis/diffbind_narrow/significant_peaks.rds")
narrow_results <- annotate_and_enrich("narrow", narrow_peaks)

# Create combined analysis summary
log_message("Creating combined analysis summary...")
combined_summary <- data.frame(
    Analysis_Type = c("Broad", "Narrow"),
    Total_Peaks = c(length(broad_peaks), length(narrow_peaks)),
    Unique_Genes = c(length(unique(broad_results$gene_list$SYMBOL)),
                    length(unique(narrow_results$gene_list$SYMBOL))),
    GO_Terms = c(nrow(broad_results$go_enrichment), 
                nrow(narrow_results$go_enrichment))
)

write.csv(combined_summary,
          "analysis/annotation_combined/combined_summary.csv",
          row.names = FALSE)

log_message("Analysis completed successfully") 