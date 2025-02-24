#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples. It identifies regions with significantly different
#' H2AK119Ub levels between conditions.
#'
#' Input files:
#' - analysis/5_peak_calling/{sample}_broad_peaks_final.broadPeak: Final peak calls for each sample
#' - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample
#'
#' Output files:
#' - analysis/7_differential_binding/: Directory containing DiffBind results
#' - analysis/7_differential_binding/: Directory containing annotated differential binding results
#'
#' Dependencies:
#' - DiffBind for differential binding analysis
#' - ChIPseeker for peak annotation
#' - TxDb.Hsapiens.UCSC.hg38.knownGene for genome annotations
#' - org.Hs.eg.db for gene ID mapping
#' - tidyverse for data manipulation
#' - rtracklayer for genomic data handling
#' - GenomicRanges for genomic interval operations
#' - clusterProfiler for enrichment analysis

# Load required libraries
suppressPackageStartupMessages({
    library(DiffBind)
    library(tidyverse)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(clusterProfiler)
    library(GenomeInfoDb)  # For renaming seqlevels
})

# Define global constants
PEAKS_DIR <- "analysis/5_peak_calling"
PEAKS_SUFFIX <- "peaks_final"
ALIGN_DIR <- "analysis/3_alignment"
OUTPUT_DIR <- commandArgs(trailingOnly = TRUE)[1] # Get output directory from command line
ANNOTATION_DIR <- file.path(OUTPUT_DIR, "../annotation_broad_deseq2") # annotation directory is one level up

# Log a message with timestamp and level
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Check if all required input files exist
validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        log_message(paste("Missing files:", paste(missing, collapse="\n")), level = "ERROR")
        stop("Missing files. Check the log for details.")
    }
    TRUE
}

# Create directories if they don't exist
create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Modify the peak file standardization:
# Since your BAM headers use names like "1", "2", … (non-UCSC style),
# remove any "chr" prefix instead of adding it.
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s", file))
        tryCatch({
            peaks <- try({
                read.table(file, header=FALSE,
                           colClasses=c("character", "numeric", "numeric", "character",
                                        "numeric", "character", "numeric", "numeric", "numeric"),
                           col.names=c("chr", "start", "end", "name", "score",
                                       "strand", "signalValue", "pValue", "qValue"))
            }, silent=TRUE)
            
            if (inherits(peaks, "try-error")) {
                peaks <- read.delim(file, header=FALSE, fill=TRUE, sep="\t",
                                    col.names=c("chr", "start", "end", "name", "score",
                                                "strand", "signalValue", "pValue", "qValue"))
            }
            
            peaks <- na.omit(peaks)
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No valid peaks found in %s after cleaning", file))
                next
            }
            
            # **Remove** any 'chr' prefix if present (so peak names match BAM files)
            peaks$chr <- gsub("^chr", "", peaks$chr)
            
            # Keep only standard chromosomes (non-UCSC style)
            standard_chroms <- c(as.character(1:22), "X", "Y", "MT")
            peaks <- peaks[peaks$chr %in% standard_chroms, ]
            
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No peaks left after filtering chromosomes in %s", file))
                next
            }
            
            write.table(peaks, file, sep="\t", quote=FALSE, 
                        row.names=FALSE, col.names=FALSE)
            
            log_message(sprintf("Processed %d peaks in %s", nrow(peaks), file))
        }, error = function(e) {
            log_message(sprintf("ERROR processing %s: %s", file, e$message), level = "ERROR")
            stop(sprintf("Failed to process peak file %s", file))
        })
    }
}

# Check and create BAM index files if missing
check_bam_indexes <- function(bam_files) {
    for(bam in bam_files) {
        if(!file.exists(paste0(bam, ".bai"))) {
            log_message(sprintf("Creating index for %s", bam), "WARNING")
            system2("samtools", args = c("index", bam))
        }
    }
}

# Perform differential binding analysis for broad peaks
perform_diffbind <- function(samples) {
    log_message("Starting differential binding analysis for broad peaks")
    
    # Create output directories
    dirs <- c(OUTPUT_DIR, ANNOTATION_DIR)
    create_dirs(dirs)
    
    # Validate input files
    validate_files(c(samples$bamReads, samples$Peaks))
    
    # Check BAM indexes
    check_bam_indexes(samples$bamReads)
    
    # Standardize chromosome names in peak files
    standardize_peak_files(samples$Peaks)
    
    # Create DiffBind object (parameters optimized for broad peaks)
    log_message("Creating DiffBind object...")
    dba_data <- dba(sampleSheet = samples,
                    minOverlap = 1,    # less stringent for broad peaks
                    peakFormat = "bed",
                    peakCaller = "broad",
                    config = data.frame(AnalysisMethod = "DBA_DESEQ2",
                                        fragmentSize = 150,
                                        doBlacklist = TRUE,
                                        RunParallel = TRUE,
                                        summits = FALSE))  # Disable summit calculation for broad peaks
    
    log_message("Initial peak counts per sample:")
    print(dba.show(dba_data))
    
    # Count reads (using parameters optimized for broad peaks)
    log_message("Counting reads...")
    dba_data <- dba.count(dba_data, 
                          bUseSummarizeOverlaps = TRUE,
                          minCount = 0,       
                          bRemoveDuplicates = TRUE,
                          score = DBA_SCORE_READS,
                          summits = FALSE)
    
    log_message("After counting - sample statistics:")
    count_info <- dba.show(dba_data)
    print(count_info)
    
    if(sum(count_info$Reads) == 0) {
        stop("No reads were counted across any samples")
    }
    
    # Normalize and perform differential analysis
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, 
                              normalize = DBA_NORM_LIB,
                              background = TRUE)
    
    dba_data <- dba.contrast(dba_data, 
                             categories = DBA_CONDITION,
                             minMembers = 2)
    
    log_message("Contrast information:")
    print(dba.show(dba_data, bContrasts = TRUE))
    
    dba_data <- dba.analyze(dba_data, 
                           method = DBA_DESEQ2)
    
    log_message("Extracting results...")
    dba_results <- dba.report(dba_data,
                             th = 1,
                             bCalled = FALSE,
                             bNormalized = TRUE)
    
    log_message("Saving results...")
    saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_analysis.rds"))
    saveRDS(dba_results, file.path(OUTPUT_DIR, "significant_peaks.rds"))
    
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, file.path(OUTPUT_DIR, "differential_peaks.csv"), row.names = FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # === Peak Annotation ===
    log_message("Annotating peaks...")
    peaks_gr <- dba_results
    
    # For annotation we need UCSC-style names. Convert the non-UCSC (e.g. "1", "2", …)
    # to UCSC style by adding the "chr" prefix.
    seqlevels(peaks_gr) <- paste0("chr", seqlevels(peaks_gr))
    
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(peaks_gr, TxDb = txdb,
                             tssRegion = c(-3000, 3000),
                             verbose = FALSE)
    
    saveRDS(peakAnno, file.path(ANNOTATION_DIR, "peak_annotation.rds"))
    write.csv(as.data.frame(peakAnno),
              file.path(ANNOTATION_DIR, "peak_annotation.csv"),
              row.names = FALSE)
    
    pdf(file.path(ANNOTATION_DIR, "annotation_plots.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    return(dba_results)
}

# Define sample list with only selected samples that have similar peak counts
# Using samples with comparable peak counts (27k-37k peaks) for more reliable comparison
samples <- data.frame(
    SampleID = c("GFP_1", "GFP_3", "YAF_2", "YAF_3"),
    Factor = rep("H2AK119Ub", 4),
    Condition = rep(c("GFP", "YAF"), each = 2),
    Replicate = c(1, 3, 2, 3),
    bamReads = file.path(ALIGN_DIR, paste0(c("GFP_1", "GFP_3", "YAF_2", "YAF_3"), ".dedup.bam")),
    Peaks = file.path(PEAKS_DIR, paste0(c("GFP_1", "GFP_3", "YAF_2", "YAF_3"), "_broad_", PEAKS_SUFFIX, ".broadPeak")),
    PeakCaller = rep("broad", 4)
)

log_message("Using selected samples with similar peak counts:")
print(samples)

# Execute the differential binding analysis
dba_results <- perform_diffbind(samples)

# Save sample sheet
sample_sheet_file <- file.path(OUTPUT_DIR, "complete_sample_sheet.csv")
write.csv(samples, sample_sheet_file, row.names = FALSE)

log_message("Successfully completed differential binding analysis")
