# Load required libraries
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
})

# Function to process peaks for a given peak type
process_peaks <- function(peak_type) {
  # Create output directories if they don't exist
  base_dir <- file.path("analysis", paste0("gene_lists_", peak_type))
  plots_dir <- file.path("analysis", paste0("plots_", peak_type))
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read the DiffBind results with error checking
  diff_peaks_file <- file.path("analysis", paste0("diffbind_", peak_type), "differential_peaks.csv")
  if (!file.exists(diff_peaks_file)) {
    stop("Input file not found: ", diff_peaks_file)
  }
  
  # Read the CSV file
  diff_peaks <- read.csv(diff_peaks_file)
  
  # Add debugging information
  message(paste0("Processing ", peak_type, " peaks"))
  message("Dimensions of input diff_peaks: ", paste(dim(diff_peaks), collapse = " x "))
  
  # Filter for YAF-enriched peaks (positive fold change and significant)
  yaf_enriched <- diff_peaks %>%
    filter(Fold > 0, FDR < 0.05)
  
  # Add debugging information
  message("Dimensions of yaf_enriched after filtering: ", paste(dim(yaf_enriched), collapse = " x "))
  message("Number of significant peaks: ", nrow(yaf_enriched))
  
  # Standardize chromosome names by adding "chr" prefix if not present
  yaf_enriched$seqnames <- ifelse(!grepl("^chr", yaf_enriched$seqnames, ignore.case = TRUE),
                            paste0("chr", yaf_enriched$seqnames),
                            yaf_enriched$seqnames)
  
  # Convert chromosome names to standard format
  yaf_enriched$seqnames <- sub("chrMT", "chrM", yaf_enriched$seqnames)
  yaf_enriched$seqnames <- sub("chr23", "chrX", yaf_enriched$seqnames)
  yaf_enriched$seqnames <- sub("chr24", "chrY", yaf_enriched$seqnames)
  
  # Add debugging information for chromosome names
  message("Unique chromosome names: ", paste(unique(yaf_enriched$seqnames), collapse = ", "))
  
  # Check if yaf_enriched has rows before proceeding
  if (nrow(yaf_enriched) == 0) {
    stop("No significant YAF-enriched peaks found")
  }
  
  # Create GRanges object for annotation
  peaks_gr <- with(yaf_enriched,
                   GRanges(seqnames = seqnames,
                          ranges = IRanges(start = start, end = end),
                          strand = "*",
                          score = Fold))
  
  # Annotate peaks
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  peakAnno <- annotatePeak(peaks_gr,
                          TxDb = txdb,
                          tssRegion = c(-3000, 3000),
                          verbose = FALSE)
  
  # Extract gene information
  genes_df <- as.data.frame(peakAnno)
  
  # Map ENTREZID to SYMBOL
  gene_symbols <- mapIds(org.Hs.eg.db,
                        keys = genes_df$geneId,
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
  
  # Create final gene list with symbols
  gene_list <- data.frame(
    ENTREZID = genes_df$geneId,
    SYMBOL = gene_symbols,
    distanceToTSS = genes_df$distanceToTSS,
    annotation = genes_df$annotation,
    fold_change = genes_df$score,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(desc(fold_change))
  
  # Save gene lists
  write.csv(gene_list,
            file = file.path(base_dir, paste0("YAF_enriched_genes_", peak_type, "_full.csv")),
            row.names = FALSE)
  
  # Save just the gene symbols
  writeLines(unique(gene_list$SYMBOL),
            file.path(base_dir, paste0("YAF_enriched_genes_", peak_type, "_symbols.txt")))
  
  # Generate plots
  pdf(file.path(plots_dir, paste0("peak_annotation_", peak_type, ".pdf")))
  print(plotAnnoBar(peakAnno))
  print(plotDistToTSS(peakAnno))
  print(upsetplot(peakAnno))
  dev.off()
  
  # Perform GO enrichment analysis
  ego <- enrichGO(gene = unique(gene_list$ENTREZID),
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  
  if (nrow(ego) > 0) {
    write.csv(as.data.frame(ego),
              file = file.path(base_dir, paste0("GO_enrichment_", peak_type, ".csv")),
              row.names = FALSE)
    
    pdf(file.path(plots_dir, paste0("GO_enrichment_", peak_type, ".pdf")))
    print(dotplot(ego, showCategory = 20))
    print(enrichMap(ego, n = 30, vertex.label.cex = 0.6))
    dev.off()
  }
  
  # Create summary
  summary_text <- c(
    paste0("YAF-enriched genes analysis summary (", peak_type, " peaks):"),
    "================================================",
    paste("Total peaks analyzed:", nrow(diff_peaks)),
    paste("YAF-enriched significant peaks:", nrow(yaf_enriched)),
    paste("Unique genes identified:", length(unique(gene_list$SYMBOL))),
    paste("GO terms enriched:", nrow(ego)),
    "\nTop 20 genes by fold change:",
    paste(head(gene_list$SYMBOL, 20), collapse = ", ")
  )
  
  writeLines(summary_text, file.path(base_dir, paste0("analysis_summary_", peak_type, ".txt")))
  
  return(list(
    gene_list = gene_list,
    peakAnno = peakAnno,
    ego = ego
  ))
}

# Process both narrow and broad peaks
narrow_results <- process_peaks("narrow")
broad_results <- process_peaks("broad")

# Create a combined summary
combined_summary <- c(
  "Combined Analysis Summary",
  "======================",
  "",
  paste("Narrow peaks unique genes:", length(unique(narrow_results$gene_list$SYMBOL))),
  paste("Broad peaks unique genes:", length(unique(broad_results$gene_list$SYMBOL))),
  "",
  "Common genes between narrow and broad peaks:",
  "-------------------------------------------",
  paste(intersect(unique(narrow_results$gene_list$SYMBOL),
                 unique(broad_results$gene_list$SYMBOL)),
        collapse = ", ")
)

writeLines(combined_summary, "analysis/gene_lists_combined/combined_analysis_summary.txt")

message("Analysis completed for both narrow and broad peaks")
