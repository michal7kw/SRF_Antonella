#!/usr/bin/env Rscript

# Load required libraries
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Set paths
h2a_peaks_dir <- "../SRF_H2AK119Ub/1_iterative_processing/analysis/peaks/"
v5_peaks_file <- "../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_peaks.narrowPeak"
diffbind_file <- "../SRF_H2AK119Ub/1_iterative_processing/analysis/diffbind/all_differential_peaks.txt"
output_dir <- "."

# Function to standardize chromosome names
standardize_chromosomes <- function(gr) {
    # Get current seqlevels
    current_levels <- seqlevels(gr)
    
    # Add 'chr' prefix if missing
    new_levels <- current_levels
    for(i in seq_along(current_levels)) {
        if(!grepl("^chr", current_levels[i]) && !grepl("^GL|^KI|^MT", current_levels[i])) {
            new_levels[i] <- paste0("chr", current_levels[i])
        }
    }
    
    # Create mapping
    names(new_levels) <- current_levels
    
    # Rename the seqlevels
    seqlevels(gr) <- new_levels
    
    # Keep only main chromosomes that exist in the data
    main_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    existing_chroms <- intersect(seqlevels(gr), main_chroms)
    
    if(length(existing_chroms) == 0) {
        stop("No standard chromosomes found in the data")
    }
    
    gr <- keepSeqlevels(gr, existing_chroms, pruning.mode="coarse")
    
    return(gr)
}

# Function to merge peaks from replicates
merge_peaks <- function(peak_files) {
    peaks_list <- lapply(peak_files, function(f) {
        message("Processing file: ", basename(f))
        peaks <- import(f)
        message("Original seqlevels: ", paste(seqlevels(peaks), collapse=", "))
        peaks <- standardize_chromosomes(peaks)
        message("Standardized seqlevels: ", paste(seqlevels(peaks), collapse=", "))
        return(peaks)
    })
    
    # Combine all peaks
    message("Combining peaks from all replicates")
    all_peaks <- unlist(GRangesList(peaks_list))
    
    # Reduce overlapping peaks
    message("Reducing overlapping peaks")
    reduced_peaks <- reduce(all_peaks)
    return(reduced_peaks)
}

# Read and standardize V5 peaks
message("Processing V5 peaks")
v5_peaks <- import(v5_peaks_file)
message("V5 original seqlevels: ", paste(seqlevels(v5_peaks), collapse=", "))
v5_peaks <- standardize_chromosomes(v5_peaks)
message("V5 standardized seqlevels: ", paste(seqlevels(v5_peaks), collapse=", "))

# Read YAF and GFP peaks separately
yaf_files <- list.files(h2a_peaks_dir, pattern="YAF_.*broadPeak$", full.names=TRUE)
gfp_files <- list.files(h2a_peaks_dir, pattern="GFP_.*broadPeak$", full.names=TRUE)

message("Processing YAF peaks")
yaf_peaks <- merge_peaks(yaf_files)
message("Processing GFP peaks")
gfp_peaks <- merge_peaks(gfp_files)

# Ensure all peak sets use the same chromosome set
common_chroms <- Reduce(intersect, list(
    seqlevels(v5_peaks),
    seqlevels(yaf_peaks),
    seqlevels(gfp_peaks)
))

message("Common chromosomes: ", paste(common_chroms, collapse=", "))

v5_peaks <- keepSeqlevels(v5_peaks, common_chroms, pruning.mode="coarse")
yaf_peaks <- keepSeqlevels(yaf_peaks, common_chroms, pruning.mode="coarse")
gfp_peaks <- keepSeqlevels(gfp_peaks, common_chroms, pruning.mode="coarse")

# Find overlaps for YAF and GFP separately
yaf_v5_overlaps <- findOverlaps(yaf_peaks, v5_peaks)
gfp_v5_overlaps <- findOverlaps(gfp_peaks, v5_peaks)

# Create overlap statistics for YAF
yaf_overlap_stats <- data.frame(
    condition = "YAF",
    h2a_peak = queryHits(yaf_v5_overlaps),
    v5_peak = subjectHits(yaf_v5_overlaps),
    v5_score = v5_peaks$score[subjectHits(yaf_v5_overlaps)]
)

# Create overlap statistics for GFP
gfp_overlap_stats <- data.frame(
    condition = "GFP",
    h2a_peak = queryHits(gfp_v5_overlaps),
    v5_peak = subjectHits(gfp_v5_overlaps),
    v5_score = v5_peaks$score[subjectHits(gfp_v5_overlaps)]
)

# Combine statistics
condition_overlap_stats <- rbind(yaf_overlap_stats, gfp_overlap_stats)

# Save condition-specific overlap statistics
write.csv(condition_overlap_stats, 
         file.path(output_dir, "condition_specific_overlaps.csv"), 
         row.names = FALSE)

# Create condition-specific overlap plots
pdf(file.path(output_dir, "condition_specific_overlaps.pdf"))
# Plot number of overlapping peaks per condition
ggplot(condition_overlap_stats, aes(x = condition)) +
    geom_bar() +
    theme_minimal() +
    labs(title = "Number of V5 overlapping peaks per condition",
         y = "Number of overlapping peaks")

# Plot V5 scores distribution by condition
ggplot(condition_overlap_stats, aes(x = condition, y = v5_score)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "V5 peak scores in overlapping regions",
         y = "V5 peak score")
dev.off()

# Annotate condition-specific overlapping regions
yaf_overlapping_peaks <- yaf_peaks[queryHits(yaf_v5_overlaps)]
gfp_overlapping_peaks <- gfp_peaks[queryHits(gfp_v5_overlaps)]

# Add condition information
mcols(yaf_overlapping_peaks)$condition <- "YAF"
mcols(gfp_overlapping_peaks)$condition <- "GFP"

# Combine peaks for annotation
all_condition_peaks <- c(yaf_overlapping_peaks, gfp_overlapping_peaks)

# Prepare TxDb object with matching chromosome naming
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) <- seqlevels(all_condition_peaks)

# Annotate peaks
condition_peak_anno <- annotatePeak(all_condition_peaks, 
                                  TxDb = txdb,
                                  tssRegion = c(-3000, 3000),
                                  annoDb = "org.Hs.eg.db")

# Save condition-specific annotation results
write.csv(as.data.frame(condition_peak_anno),
         file.path(output_dir, "condition_specific_peak_annotation.csv"),
         row.names = FALSE)

# Create condition-specific annotation visualization
pdf(file.path(output_dir, "condition_specific_annotation_plots.pdf"))
plotAnnoPie(condition_peak_anno)
plotDistToTSS(condition_peak_anno, facet = "condition")
dev.off()

# Generate condition-specific summary statistics
cat("Condition-specific Summary Statistics:\n",
    "Total YAF peaks:", length(yaf_peaks), "\n",
    "Total GFP peaks:", length(gfp_peaks), "\n",
    "Total V5 peaks:", length(v5_peaks), "\n",
    "YAF regions overlapping with V5:", length(unique(queryHits(yaf_v5_overlaps))), "\n",
    "GFP regions overlapping with V5:", length(unique(queryHits(gfp_v5_overlaps))), "\n",
    file = file.path(output_dir, "condition_specific_summary.txt"))

# Export overlapping regions as BED files for visualization
export(yaf_peaks[queryHits(yaf_v5_overlaps)], 
       file.path(output_dir, "yaf_v5_overlapping.bed"))
export(gfp_peaks[queryHits(gfp_v5_overlaps)], 
       file.path(output_dir, "gfp_v5_overlapping.bed"))

# Export all peaks for reference
export(yaf_peaks, file.path(output_dir, "yaf_all_peaks.bed"))
export(gfp_peaks, file.path(output_dir, "gfp_all_peaks.bed"))
export(v5_peaks, file.path(output_dir, "v5_all_peaks.bed"))
