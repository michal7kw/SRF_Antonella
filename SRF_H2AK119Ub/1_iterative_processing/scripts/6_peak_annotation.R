# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# Function to standardize chromosome names
standardize_chromosomes <- function(gr) {
    # First, check if chromosomes have 'chr' prefix
    current_chroms <- seqlevels(gr)
    
    # Add 'chr' prefix if missing
    if (!all(grepl("^chr", current_chroms))) {
        seqlevels(gr) <- paste0("chr", seqlevels(gr))
    }
    
    # Define standard chromosomes
    standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    
    # Keep only standard chromosomes
    gr <- gr[seqnames(gr) %in% standard_chroms]
    
    # Set standard chromosome levels
    seqlevels(gr) <- standard_chroms[standard_chroms %in% seqlevels(gr)]
    
    # Ensure genome info is set
    genome(gr) <- "hg38"
    
    return(gr)
}

# Function to perform annotation analysis
perform_annotation <- function(peak_type) {
    # Create output directories
    base_dir <- file.path("analysis", paste0("annotation_", peak_type))
    dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    # Load differential binding results
    print(paste("Loading", peak_type, "peak differential binding results..."))
    diffbind_dir <- file.path("analysis", paste0("diffbind_", peak_type))
    all_peaks <- readRDS(file.path(diffbind_dir, "all_peaks.rds"))
    sig_peaks <- readRDS(file.path(diffbind_dir, "significant_peaks.rds"))
    
    # Standardize chromosome names
    print("Standardizing chromosome names...")
    all_peaks <- standardize_chromosomes(all_peaks)
    sig_peaks <- standardize_chromosomes(sig_peaks)
    
    # Get TxDb object
    print("Loading TxDb...")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    # Annotate peaks
    print("Annotating peaks...")
    peak_anno <- annotatePeak(sig_peaks, 
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            verbose = FALSE)
    
    # Save annotation results
    print("Saving annotation results...")
    saveRDS(peak_anno, file.path(base_dir, "peak_annotation.rds"))
    write.csv(as.data.frame(peak_anno), 
             file.path(base_dir, "tables", "peak_annotation.csv"),
             row.names = FALSE)
    
    # Generate annotation plots
    print("Generating annotation plots...")
    pdf(file.path(base_dir, "figures", "annotation_plots.pdf"))
    print(plotAnnoBar(peak_anno))
    print(plotDistToTSS(peak_anno))
    print(upsetplot(peak_anno))
    dev.off()
    
    # Perform GO enrichment analysis
    print("Performing GO enrichment analysis...")
    genes <- unique(peak_anno@anno$geneId)
    ego <- tryCatch({
        enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
    }, error = function(e) {
        print(paste("Error in GO enrichment:", e))
        return(NULL)
    })
    
    if (!is.null(ego) && nrow(ego) > 0) {
        # Save GO results
        saveRDS(ego, file.path(base_dir, "go_enrichment.rds"))
        write.csv(as.data.frame(ego),
                 file.path(base_dir, "tables", "go_enrichment.csv"),
                 row.names = FALSE)
        
        # Generate GO plots
        pdf(file.path(base_dir, "figures", "go_enrichment_plots.pdf"))
        print(dotplot(ego, showCategory = 20))
        print(enrichMap(ego, n = 30, vertex.label.cex = 0.6))
        dev.off()
    }
    
    # Perform KEGG pathway analysis
    print("Performing KEGG pathway analysis...")
    entrez_ids <- mapIds(org.Hs.eg.db,
                        keys = genes,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
    
    if (!is.null(entrez_ids) && length(entrez_ids) > 0) {
        ekegg <- enrichKEGG(gene = unique(na.omit(entrez_ids)),
                           organism = "hsa",
                           pvalueCutoff = 0.05)
        
        if (nrow(ekegg) > 0) {
            saveRDS(ekegg, file.path(base_dir, "kegg_enrichment.rds"))
            write.csv(as.data.frame(ekegg),
                     file.path(base_dir, "tables", "kegg_enrichment.csv"),
                     row.names = FALSE)
            
            pdf(file.path(base_dir, "figures", "kegg_enrichment_plots.pdf"))
            print(dotplot(ekegg, showCategory = 20))
            dev.off()
        }
    }
    
    print(paste("Completed annotation analysis for", peak_type, "peaks"))
}

# Run annotation for both narrow and broad peaks
perform_annotation("narrow")
perform_annotation("broad")

print("Peak annotation analysis completed for both narrow and broad peaks")
