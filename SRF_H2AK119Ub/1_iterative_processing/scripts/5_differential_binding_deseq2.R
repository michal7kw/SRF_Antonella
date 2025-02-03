#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples using DESeq2. It identifies genomic regions with
#' significantly different H2AK119Ub levels between conditions.
#'
#' Input files:
#' - BAM files in analysis/aligned/
#'   - YAF_[1-3].dedup.bam: Aligned reads from YAF samples
#'   - GFP_[1-3].dedup.bam: Aligned reads from GFP samples
#' - Peak files in analysis/peaks2_improved/
#'   - YAF_[1-3]_broad_peaks_final.broadPeak: Called peaks from YAF samples
#'   - GFP_[1-3]_broad_peaks_final.broadPeak: Called peaks from GFP samples
#'
#' Output files in analysis/diffbind_broad_deseq2/:
#' - consensus_peaks.bed: Merged peaks from all samples used for counting
#' - differential_peaks_deseq2.csv: Table of differential peaks with statistics
#' - deseq2_analysis.rds: Complete DESeq2 analysis object
#' - MA_plot.pdf: MA plot showing fold changes vs mean normalized counts
#'
#' Output files in analysis/annotation_broad_deseq2/:
#' - significant_peaks_annotation.csv: Annotation of significant peaks
#' - annotation_plots.pdf: Visualizations of peak annotations
#'
#' Dependencies:
#' - GenomicRanges and GenomicAlignments for genomic data handling
#' - DESeq2 for differential analysis
#' - ChIPseeker for peak annotation
#' - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
#' - org.Hs.eg.db for gene ID mapping
#' - tidyverse for data manipulation
#' - rtracklayer for genomic data import/export

# Load required libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomicAlignments)
    library(rtracklayer)
    library(DESeq2)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(tidyverse)
})

# Define global constants for input/output directories
PEAKS_DIR <- "analysis/peaks2_improved"
OUTPUT_DIR <- "analysis/diffbind_broad_deseq2"
ANNOTATION_DIR <- "analysis/annotation_broad_deseq2"

# Utility function for logging messages with timestamp
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOTATION_DIR, recursive = TRUE, showWarnings = FALSE)

# Define sample metadata including file paths
samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Condition = rep(c("YAF", "GFP"), each=3),
    Replicate = rep(1:3, 2),
    PeakFile = file.path(PEAKS_DIR, paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)), 
                                          "_broad_peaks_final.broadPeak")),
    BamFile = file.path("analysis/aligned", paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)), 
                                                  ".dedup.bam"))
)

# Read and merge peaks from all samples with diagnostic output
log_message("Reading and merging peaks...")
peak_list <- list()
for (i in seq_along(samples$PeakFile)) {
    f <- samples$PeakFile[i]
    log_message(sprintf("Reading peaks from %s", basename(f)))
    peaks <- read.table(f, col.names=c("chr", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue"))
    log_message(sprintf("Found %d peaks in %s", nrow(peaks), basename(f)))
    peak_list[[i]] <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
}

# Create consensus peaks by merging nearby peaks (within 1kb)
log_message("Creating consensus peaks...")
all_peaks <- unlist(GRangesList(peak_list))
log_message(sprintf("Total peaks before merging: %d", length(all_peaks)))

consensus_peaks <- GenomicRanges::reduce(all_peaks, min.gapwidth=1000)
log_message(sprintf("Created consensus peak set with %d regions", length(consensus_peaks)))

# Save consensus peaks in BED format for visualization
log_message("Saving consensus peaks...")
consensus_bed <- data.frame(chr=seqnames(consensus_peaks),
                           start=start(consensus_peaks)-1,
                           end=end(consensus_peaks))
write.table(consensus_bed,
            file=file.path(OUTPUT_DIR, "consensus_peaks.bed"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Count reads in consensus peaks with strict parameters
log_message("Counting reads in consensus peaks...")

# Verify BAM index files exist
for (bam in samples$BamFile) {
    log_message(sprintf("Checking BAM file: %s", basename(bam)))
    if (!file.exists(paste0(bam, ".bai"))) {
        stop(sprintf("BAM index not found for %s", bam))
    }
}

# Set up counting parameters for proper paired-end read handling
bam_files <- BamFileList(samples$BamFile)
count_params <- ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                             isUnmappedQuery=FALSE,
                                             isNotPassingQualityControls=FALSE))

# Count reads using strict intersection mode
se <- summarizeOverlaps(features=consensus_peaks,
                       reads=bam_files,
                       mode="IntersectionStrict",
                       singleEnd=FALSE,
                       ignore.strand=TRUE,
                       inter.feature=TRUE,
                       param=count_params)

# Verify read counting worked
log_message("Checking read counts...")
total_counts <- colSums(assay(se))
for (i in seq_along(total_counts)) {
    log_message(sprintf("%s: %d reads counted", samples$SampleID[i], total_counts[i]))
}

if (all(total_counts == 0)) {
    stop("No reads counted in any sample. Check BAM files and consensus peaks.")
}

# Set up and run DESeq2 analysis
log_message("Creating DESeq2 object...")
colData(se)$condition <- factor(samples$Condition)
colData(se)$replicate <- factor(samples$Replicate)
dds <- DESeqDataSet(se, design = ~ condition)

# Filter low count regions (at least 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
log_message(sprintf("Retained %d regions after filtering", nrow(dds)))

# Run DESeq2 differential analysis
log_message("Running DESeq2 analysis...")
dds <- DESeq(dds)

# Extract and order results by adjusted p-value
log_message("Extracting results...")
res <- results(dds, contrast=c("condition", "YAF", "GFP"))
res <- res[order(res$padj),]

# Add genomic coordinates to results
peak_coords <- as.data.frame(consensus_peaks[keep])[,1:3]
res_df <- as.data.frame(res)
res_df <- cbind(peak_coords, res_df)

# Save all analysis results
log_message("Saving results...")
write.csv(res_df, file.path(OUTPUT_DIR, "differential_peaks_deseq2.csv"))
saveRDS(dds, file.path(OUTPUT_DIR, "deseq2_analysis.rds"))

# Generate MA plot
pdf(file.path(OUTPUT_DIR, "MA_plot.pdf"))
plotMA(res, ylim=c(-5,5))
dev.off()

# Annotate significant peaks (FDR < 0.05 and |log2FC| > 1)
log_message("Annotating significant peaks...")
sig_peaks <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
if(nrow(sig_peaks) > 0) {
    sig_gr <- makeGRangesFromDataFrame(sig_peaks)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(sig_gr, TxDb=txdb,
                            tssRegion=c(-3000, 3000),
                            verbose=FALSE)
    
    # Save annotation results and generate visualization plots
    write.csv(as.data.frame(peakAnno), 
              file.path(ANNOTATION_DIR, "significant_peaks_annotation.csv"),
              row.names=FALSE)
    
    pdf(file.path(ANNOTATION_DIR, "annotation_plots.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
}

# Print final summary statistics
log_message("Analysis complete!")
log_message(sprintf("Total peaks analyzed: %d", nrow(res_df)))
log_message(sprintf("Significant peaks (FDR < 0.05, |log2FC| > 1): %d", nrow(sig_peaks)))
