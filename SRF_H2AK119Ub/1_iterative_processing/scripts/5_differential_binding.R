# Load required libraries
library(DiffBind)
library(tidyverse)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(clusterProfiler)

# Create necessary directories
dir.create("analysis/diffbind_narrow", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/annotation_narrow", recursive = TRUE, showWarnings = FALSE)

# Create sample sheet for DiffBind
samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Factor = rep("H2AK119Ub", 6),
    Condition = rep(c("YAF", "GFP"), each=3),
    Replicate = rep(1:3, 2),
    bamReads = file.path("analysis/aligned", 
                        c(paste0("YAF_", 1:3, ".dedup.bam"),
                          paste0("GFP_", 1:3, ".dedup.bam"))),
    Peaks = file.path("analysis/peaks",
                     c(paste0("YAF_", 1:3, "_narrow_peaks.narrowPeak"),
                       paste0("GFP_", 1:3, "_narrow_peaks.narrowPeak"))),
    PeakCaller = rep("narrow", 6)
)

# Add error checking for file existence
check_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
}

# Check if all files exist
print("Checking input files...")
check_files(samples$bamReads)
check_files(samples$Peaks)

# Create DiffBind object
print("Creating DiffBind object...")
dba_data <- dba(sampleSheet=samples,
                minOverlap=2,
                peakFormat="bed",
                peakCaller="narrow",
                config=data.frame(AnalysisMethod="max",
                                fragmentSize=150,
                                doBlacklist=TRUE))

# Count reads
print("Counting reads...")
dba_data <- dba.count(dba_data)

# Normalize and perform differential analysis
print("Performing differential analysis...")
dba_data <- dba.normalize(dba_data)
dba_data <- dba.contrast(dba_data, 
                        categories=DBA_CONDITION,
                        minMembers=2)
dba_data <- dba.analyze(dba_data)

# Extract results
print("Extracting results...")
dba_results <- dba.report(dba_data)

# Save complete analysis
print("Saving results...")
saveRDS(dba_data, "analysis/diffbind_narrow/all_peaks.rds")
saveRDS(dba_results, "analysis/diffbind_narrow/significant_peaks.rds")

# Convert to data frame and save
df_results <- as.data.frame(dba_results)
write.csv(df_results, "analysis/diffbind_narrow/differential_peaks.csv", row.names=FALSE)

# Annotate peaks
print("Annotating peaks...")
peaks_gr <- dba_results
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(peaks_gr, TxDb=txdb,
                        tssRegion=c(-3000, 3000),
                        verbose=FALSE)

# Save annotation results
print("Saving annotation results...")
saveRDS(peakAnno, "analysis/annotation_narrow/peak_annotation.rds")
write.csv(as.data.frame(peakAnno), 
         "analysis/annotation_narrow/peak_annotation.csv",
         row.names=FALSE)

# Generate annotation plots
pdf("analysis/annotation_narrow/annotation_plots.pdf")
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)
dev.off()

# Perform GO enrichment analysis
genes <- unique(peakAnno@anno$geneId)
ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

# Save enrichment results
if (nrow(ego) > 0) {
    saveRDS(ego, "analysis/annotation_narrow/go_enrichment.rds")
    write.csv(as.data.frame(ego), 
             "analysis/annotation_narrow/go_enrichment.csv",
             row.names=FALSE)
    
    # Generate enrichment plot
    pdf("analysis/annotation_narrow/go_enrichment_plot.pdf")
    dotplot(ego, showCategory=20)
    dev.off()
}

print("Differential binding analysis for narrow peaks completed")