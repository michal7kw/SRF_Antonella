#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples. It identifies regions with significantly different
#' H2AK119Ub levels between conditions.
#'
#' Input files:
#' - BAM files in analysis/aligned/ (or analysis/aligned if that is intended)
#'   - YAF_[1-3].dedup.bam: Aligned reads from YAF samples
#'   - GFP_[1-3].dedup.bam: Aligned reads from GFP samples
#' - Peak files in analysis/peaks/ (or analysis/peaks if that is intended)
#'   - YAF_[1-3]_broad_peaks_final.broadPeak: Called peaks from YAF samples
#'   - GFP_[1-3]_broad_peaks_final.broadPeak: Called peaks from GFP samples
#'
#' Output files in analysis/diffbind_broad/:
#' - dba_analysis.rds: Complete DiffBind analysis object
#' - significant_peaks.rds: GRanges object with differential peaks
#' - differential_peaks.csv: Table of differential peaks with statistics
#' - complete_sample_sheet.csv: Sample metadata used in analysis
#'
#' Output files in analysis/annotation_broad/:
#' - peak_annotation.rds: ChIPseeker annotation object
#' - peak_annotation.csv: Table with peak annotations
#' - annotation_plots.pdf: Visualizations of peak annotations
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
# <<< Adjust these paths if needed >>>
PEAKS_DIR <- "analysis/peaks"  # change to "analysis/peaks" if that is intended
PEAKS_SUFFIX <- "peaks_final"
OUTPUT_DIR <- "analysis/diffbind_broad"
ANNOTATION_DIR <- "analysis/annotation_broad"

# Log a message with timestamp and level
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Check if all required input files exist
validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
    TRUE
}

# Create directories if they don't exist
create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Standardize peak files so that their chromosome names match the BAM files.
# Since the BAM headers use names like "1", "2", … (non-UCSC style), remove any "chr" prefix.
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s", file))
        tryCatch({
            peaks <- try({
                read.table(file, header = FALSE,
                           colClasses = c("character", "numeric", "numeric", "character",
                                          "numeric", "character", "numeric", "numeric", "numeric"),
                           col.names = c("chr", "start", "end", "name", "score",
                                         "strand", "signalValue", "pValue", "qValue"))
            }, silent = TRUE)
            
            if (inherits(peaks, "try-error")) {
                peaks <- read.delim(file, header = FALSE, fill = TRUE, sep = "\t",
                                    col.names = c("chr", "start", "end", "name", "score",
                                                  "strand", "signalValue", "pValue", "qValue"))
            }
            
            peaks <- na.omit(peaks)
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No valid peaks found in %s after cleaning", file))
                next
            }
            
            # Remove any 'chr' prefix so that peaks match the BAM file chromosome names.
            peaks$chr <- gsub("^chr", "", peaks$chr)
            
            # Keep only standard chromosomes (non-UCSC style)
            standard_chroms <- c(as.character(1:22), "X", "Y", "MT")
            peaks <- peaks[peaks$chr %in% standard_chroms, ]
            
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No peaks left after filtering chromosomes in %s", file))
                next
            }
            
            write.table(peaks, file, sep = "\t", quote = FALSE, 
                        row.names = FALSE, col.names = FALSE)
            log_message(sprintf("Processed %d peaks in %s", nrow(peaks), file))
        }, error = function(e) {
            log_message(sprintf("ERROR processing %s: %s", file, e$message))
            stop(sprintf("Failed to process peak file %s", file))
        })
    }
}

# Check and create BAM index files if missing
check_bam_indexes <- function(bam_files) {
    for (bam in bam_files) {
        if (!file.exists(paste0(bam, ".bai"))) {
            log_message(sprintf("Creating index for %s", bam), "WARNING")
            system2("samtools", args = c("index", bam))
        }
    }
}

# Perform differential binding analysis for broad peaks
perform_diffbind <- function() {
    log_message("Starting differential binding analysis for broad peaks")
    
    # Create output directories
    dirs <- c(OUTPUT_DIR, ANNOTATION_DIR, "logs")
    create_dirs(dirs)
    
    # Define a sample sheet for the selected samples.
    # (If you want to analyze all samples, adjust this data frame accordingly.)
    samples <- data.frame(
        SampleID = c("YAF_2", "YAF_3", "GFP_1", "GFP_3"),
        Factor = rep("H2AK119Ub", 4),
        Condition = rep(c("YAF", "GFP"), each = 2),
        Replicate = c(2:3, 1, 3),
        bamReads = file.path("analysis/aligned", 
                             paste0(c("YAF_2", "YAF_3", "GFP_1", "GFP_3"), ".dedup.bam")),
        Peaks = file.path(PEAKS_DIR,
                          paste0(c("YAF_2", "YAF_3", "GFP_1", "GFP_3"), "_broad_", PEAKS_SUFFIX, ".broadPeak")),
        PeakCaller = rep("broad", 4)
    )
    
    log_message("Using selected samples with similar peak counts:")
    log_message("YAF_2: ~35k peaks")
    log_message("YAF_3: ~37k peaks")
    log_message("GFP_1: ~29k peaks")
    log_message("GFP_3: ~27k peaks")
    
    # Validate input files
    validate_files(c(samples$bamReads, samples$Peaks))
    
    # Check BAM indexes
    check_bam_indexes(samples$bamReads)
    
    # Standardize chromosome names in peak files (once)
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
    
    if (sum(count_info$Reads) == 0) {
        stop("No reads were counted across any samples")
    }
    
    # Normalize and perform differential analysis
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, 
                              normalize = DBA_NORM_LIB,
                              background = TRUE)
    
    dba_data <- dba.contrast(dba_data, 
                             categories = DBA_CONDITION,
                             minMembers = 2)  # require at least 2 replicates per condition
    
    log_message("Contrast information:")
    print(dba.show(dba_data, bContrasts = TRUE))
    
    dba_data <- dba.analyze(dba_data, 
                           method = DBA_DESEQ2,
                           bFullLibrarySize = TRUE,
                           bTagwise = FALSE)
    
    log_message("Extracting results...")
    dba_results <- dba.report(dba_data,
                              th = 1,   # no threshold filtering
                              bCalled = FALSE,
                              bNormalized = TRUE,
                              bSummarized = TRUE,
                              bControlSubtracted = TRUE)
    
    log_message("Saving results...")
    saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_analysis.rds"))
    saveRDS(dba_results, file.path(OUTPUT_DIR, "significant_peaks.rds"))
    
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, file.path(OUTPUT_DIR, "differential_peaks.csv"), row.names = FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # === Peak Annotation ===
    log_message("Annotating peaks...")
    peaks_gr <- dba_results
    
    # For annotation we need UCSC-style names.
    # Convert the non-UCSC style (e.g., "1", "2", …) to UCSC by adding "chr".
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

# === Main Execution ===

# Option 1: If you wish to run the analysis only on the selected samples,
# then simply run the following. (The sample sheet is defined inside perform_diffbind.)
dba_results <- perform_diffbind()

# Option 2: If you wish to have a complete sample sheet for all samples (e.g., for record keeping),
# define it once and save it. (Make sure it uses the same peak standardization as above.)
complete_samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Factor = rep("H2AK119Ub", 6),
    Condition = rep(c("YAF", "GFP"), each = 3),
    Replicate = rep(1:3, 2),
    bamReads = file.path("analysis/aligned", 
                         paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)), ".dedup.bam")),
    Peaks = file.path(PEAKS_DIR,
                      paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)), "_broad_", PEAKS_SUFFIX, ".broadPeak")),
    PeakCaller = rep("broad", 6)
)

# Validate and standardize peaks for the complete sample sheet if desired
validate_files(c(complete_samples$bamReads, complete_samples$Peaks))
standardize_peak_files(complete_samples$Peaks)
write.csv(complete_samples, file.path(OUTPUT_DIR, "complete_sample_sheet.csv"), row.names = FALSE)

log_message("Successfully completed differential binding analysis")
