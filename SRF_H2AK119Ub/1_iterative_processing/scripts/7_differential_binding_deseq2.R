#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples using DESeq2. It identifies genomic regions with
#' significantly different H2AK119Ub levels between conditions.
#'
#' Input files:
#' - analysis/5_peak_calling/{sample}_broad_peaks_final.broadPeak: Final peak calls for each sample
#' - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample
#'
#' Output files:
#' - analysis/7_differential_binding_deseq2/: Directory containing DESeq2 results and annotated results
#'   - consensus_peaks.bed: Merged peaks from all samples used for counting
#'   - differential_peaks_deseq2.csv: Table of differential peaks with statistics
#'   - deseq2_analysis.rds: Complete DESeq2 analysis object
#'   - MA_plot.pdf: MA plot showing fold changes vs mean normalized counts
#' - analysis/annotation_broad_deseq2/:
#'   - significant_peaks_annotation.csv: Annotation of significant peaks
#'   - annotation_plots.pdf: Visualizations of peak annotations
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

# Define command line argument for output directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Output directory argument is missing. Please provide the output directory path.")
}
OUTPUT_DIR <- args[1]

# Define global constants for input directories and annotation directory
PEAKS_DIR <- "analysis/5_peak_calling"
ALIGN_DIR <- "analysis/3_alignment"
ANNOTATION_DIR <- file.path(OUTPUT_DIR, "../annotation_broad_deseq2") # annotation directory is one level up

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
    SampleID = c("GFP_1", "GFP_3", "YAF_2", "YAF_3"),
    Condition = rep(c("GFP", "YAF"), each = 2),
    Replicate = c(1, 3, 2, 3),
    PeakFile = file.path(PEAKS_DIR, paste0(c("GFP_1", "GFP_3", "YAF_2", "YAF_3"), 
                                          "_broad_peaks_final.broadPeak")),
    BamFile = file.path(ALIGN_DIR, paste0(c("GFP_1", "GFP_3", "YAF_2", "YAF_3"), 
                                                  ".dedup.bam"))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ AND MERGE PEAKS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_message("Reading and merging peaks...")

# Use a list to store GRanges objects from each peak file.
peak_list <- list()

# Extract Seqinfo from one BAM file for consistency.
bam_si <- seqinfo(BamFile(samples$BamFile[1]))

for (i in seq_along(samples$PeakFile)) {
    f <- samples$PeakFile[i]
    log_message(sprintf("Reading peaks from %s", basename(f)))
    # Read the peak file (assuming no header and 9 columns)
    peaks <- read.table(f, col.names = c("chr", "start", "end", "name", "score",
                                           "strand", "signalValue", "pValue", "qValue"))
    log_message(sprintf("Found %d peaks in %s", nrow(peaks), basename(f)))
    
    # Remove any "chr" prefix so that peak seqlevels match BAM files
    peaks$chr <- gsub("^chr", "", peaks$chr)
    
    # Convert to GRanges
    gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
    
    # Restrict to sequence levels that appear in the BAM file
    common_levels <- intersect(seqlevels(gr), seqlevels(bam_si))
    gr <- keepSeqlevels(gr, common_levels, pruning.mode = "coarse")
    # Set the Seqinfo to that from the BAM file (so all GRanges objects have common seqinfo)
    seqinfo(gr) <- bam_si[common_levels]
    
    peak_list[[i]] <- gr
}

# Merge all peaks into a single GRanges object using a loop
all_peaks <- peak_list[[1]]
if (length(peak_list) > 1) {
    for (i in 2:length(peak_list)) {
        all_peaks <- c(all_peaks, peak_list[[i]])
    }
}
log_message(sprintf("Total peaks before merging: %d", length(all_peaks)))

# Create consensus peaks by merging nearby peaks (within 1kb)
# First convert to ranges without metadata to avoid reduce() issues
ranges_only <- GRanges(seqnames(all_peaks), 
                      ranges(all_peaks),
                      seqinfo = seqinfo(all_peaks))
consensus_peaks <- GenomicRanges::reduce(ranges_only, min.gapwidth = 1000)
log_message(sprintf("Created consensus peak set with %d regions", length(consensus_peaks)))

# Save consensus peaks in BED format for visualization
log_message("Saving consensus peaks...")
consensus_bed <- data.frame(chr = as.character(seqnames(consensus_peaks)),
                            start = start(consensus_peaks) - 1,
                            end = end(consensus_peaks))
write.table(consensus_bed,
            file = file.path(OUTPUT_DIR, "consensus_peaks.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ COUNTING WITH SUMMARIZEOverlaps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
count_params <- ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                                  isUnmappedQuery = FALSE,
                                                  isNotPassingQualityControls = FALSE))

# Count reads using strict intersection mode
se <- summarizeOverlaps(features = consensus_peaks,
                        reads = bam_files,
                        mode = "IntersectionStrict",
                        singleEnd = FALSE,
                        ignore.strand = TRUE,
                        inter.feature = TRUE,
                        param = count_params)

# Verify read counting worked
log_message("Checking read counts...")
total_counts <- colSums(assay(se))
for (i in seq_along(total_counts)) {
    log_message(sprintf("%s: %d reads counted", samples$SampleID[i], total_counts[i]))
}

if (all(total_counts == 0)) {
    stop("No reads counted in any sample. Check BAM files and consensus peaks.")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESeq2 ANALYSIS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_message("Creating DESeq2 object...")
colData(se)$condition <- factor(samples$Condition)
colData(se)$replicate <- factor(samples$Replicate)
dds <- DESeqDataSet(se, design = ~ condition)

# Filter low count regions (at least 10 reads total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
log_message(sprintf("Retained %d regions after filtering", nrow(dds)))

# Run DESeq2 differential analysis
log_message("Running DESeq2 analysis...")
dds <- DESeq(dds)

# Extract and order results by adjusted p-value
log_message("Extracting results...")
res <- results(dds, contrast = c("condition", "YAF", "GFP"))
res <- res[order(res$padj), ]

# Add genomic coordinates to results.
# Extract the consensus peaks corresponding to the rows that passed filtering.
peak_coords <- as.data.frame(consensus_peaks[keep])[, 1:3]
res_df <- as.data.frame(res)
res_df <- cbind(peak_coords, res_df)

# Save all analysis results
log_message("Saving results...")
write.csv(res_df, file = file.path(OUTPUT_DIR, "differential_peaks_deseq2.csv"), row.names = FALSE)
saveRDS(dds, file = file.path(OUTPUT_DIR, "deseq2_analysis.rds"))

# Generate MA plot
pdf(file = file.path(OUTPUT_DIR, "MA_plot.pdf"))
plotMA(res, ylim = c(-5, 5))
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PEAK ANNOTATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log_message("Annotating significant peaks...")
sig_peaks <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
if (nrow(sig_peaks) > 0) {
    # Create GRanges from the significant peaks.
    # Ensure column names match expected fields
    names(sig_peaks)[names(sig_peaks) == "chr"] <- "seqnames"  # Rename chr to seqnames
    sig_gr <- makeGRangesFromDataFrame(sig_peaks,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE)
    
    # Convert to UCSC-style by adding "chr" if needed.
    seqlevels(sig_gr) <- paste0("chr", seqlevels(sig_gr))
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(sig_gr, TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            verbose = FALSE)
    
    write.csv(as.data.frame(peakAnno), 
              file = file.path(ANNOTATION_DIR, "significant_peaks_annotation.csv"),
              row.names = FALSE)
    
    pdf(file = file.path(ANNOTATION_DIR, "annotation_plots.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
}

# Final summary statistics
log_message("Analysis complete!")
log_message(sprintf("Total peaks analyzed: %d", nrow(res_df)))
log_message(sprintf("Significant peaks (FDR < 0.05, |log2FC| > 1): %d", nrow(sig_peaks)))
