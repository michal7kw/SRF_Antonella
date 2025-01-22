# Load required libraries
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
})

# Create output directories if they don't exist
dir.create("analysis/gene_lists", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/plots", recursive = TRUE, showWarnings = FALSE)

# Read the DiffBind results with error checking
diff_peaks_file <- "analysis/diffbind/all_differential_peaks.txt"
if (!file.exists(diff_peaks_file)) {
  stop("Input file not found: ", diff_peaks_file)
}

# Read the tab-delimited file instead of RDS
diff_peaks <- read.table(diff_peaks_file, header = TRUE, sep = "\t")

# Add debugging information
message("Dimensions of input diff_peaks: ", paste(dim(diff_peaks), collapse = " x "))

# Filter for YAF-enriched peaks (positive fold change and significant)
yaf_enriched <- diff_peaks %>%
  filter(Fold > 0, FDR < 0.05)

# Add debugging information
message("Dimensions of yaf_enriched after filtering: ", paste(dim(yaf_enriched), collapse = " x "))
message("Number of significant peaks: ", nrow(yaf_enriched))

# Standardize chromosome names by adding "chr" prefix if not present
yaf_enriched$Chr <- ifelse(!grepl("^chr", yaf_enriched$Chr, ignore.case = TRUE),
                          paste0("chr", yaf_enriched$Chr),
                          yaf_enriched$Chr)

# Convert chromosome names to standard format
yaf_enriched$Chr <- sub("chrMT", "chrM", yaf_enriched$Chr)
yaf_enriched$Chr <- sub("chr23", "chrX", yaf_enriched$Chr)
yaf_enriched$Chr <- sub("chr24", "chrY", yaf_enriched$Chr)

# Add debugging information for chromosome names
message("Unique chromosome names: ", paste(unique(yaf_enriched$Chr), collapse = ", "))

# Check if yaf_enriched has rows before proceeding
if (nrow(yaf_enriched) == 0) {
  stop("No significant YAF-enriched peaks found")
}

# Convert to GRanges object with standardized chromosome names
yaf_peaks_gr <- GRanges(
  seqnames = yaf_enriched$Chr,
  ranges = IRanges(start = yaf_enriched$Start, end = yaf_enriched$End),
  strand = "*",
  fold_change = yaf_enriched$Fold,
  FDR = yaf_enriched$FDR
)

# Add debugging information for GRanges object
message("GRanges seqnames: ", paste(unique(seqnames(yaf_peaks_gr)), collapse = ", "))

# Annotate peaks with nearby genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(yaf_peaks_gr, 
                        TxDb = txdb,
                        tssRegion = c(-3000, 3000),
                        annoDb = "org.Hs.eg.db")

# Convert annotation to data frame
anno_df <- as.data.frame(peakAnno)

# Add gene symbols with error handling
anno_df$gene_symbol <- tryCatch({
  mapIds(org.Hs.eg.db,
         keys = anno_df$geneId,
         column = "SYMBOL",
         keytype = "ENTREZID",
         multiVals = "first")
}, error = function(e) {
  message("Error mapping gene symbols: ", e$message)
  return(NA)
})

# Create simplified gene list with key information
simplified_list <- anno_df %>%
  select(gene_symbol, annotation, fold_change = fold_change, FDR,
         distance_to_TSS = distanceToTSS) %>%
  arrange(desc(fold_change)) %>%
  filter(!is.na(gene_symbol))  # Remove entries without gene symbols

# Save results
write.table(anno_df,
            "analysis/gene_lists/YAF_enriched_H2AK119Ub_full_annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(simplified_list,
            "analysis/gene_lists/YAF_enriched_H2AK119Ub_genes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Perform GO enrichment analysis
gene_list <- unique(na.omit(anno_df$geneId))
if (length(gene_list) > 0) {
  tryCatch({
    go_bp <- enrichGO(gene = gene_list,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
    
    if (nrow(go_bp) > 0) {
      # Save GO results
      write.table(as.data.frame(go_bp),
                  "analysis/gene_lists/YAF_enriched_GO_terms.txt",
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      # Create GO dotplot
      pdf("analysis/plots/GO_enrichment_dotplot.pdf", width = 10, height = 8)
      print(dotplot(go_bp, showCategory = 20))
      dev.off()
    } else {
      message("No significant GO terms found")
    }
  }, error = function(e) {
    message("Error in GO enrichment analysis: ", e$message)
  })
}

# Print summary statistics
summary_text <- sprintf(
  "Summary of YAF-enriched H2AK119Ub regions:\n
  Total significant peaks: %d\n
  Total annotated genes: %d\n
  Number of significant GO terms: %d\n
  Analysis completed: %s",
  nrow(yaf_enriched),
  length(unique(na.omit(anno_df$gene_symbol))),
  if(exists("go_bp")) nrow(go_bp) else 0,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(summary_text, "analysis/gene_lists/analysis_summary.txt")
