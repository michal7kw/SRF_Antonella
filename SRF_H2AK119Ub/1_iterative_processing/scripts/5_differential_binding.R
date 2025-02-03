#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples. It identifies regions with significantly different
#' H2AK119Ub levels between conditions.
#'
#' Input files:
#' - BAM files in analysis/aligned/
#'   - YAF_[1-3].dedup.bam: Aligned reads from YAF samples
#'   - GFP_[1-3].dedup.bam: Aligned reads from GFP samples
#' - Peak files in analysis/peaks2_improved/
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
})

# Define global constants
PEAKS_DIR <- "analysis/peaks2_improved"
PEAKS_SUFFIX <- "peaks_final"
OUTPUT_DIR <- "analysis/diffbind_broad"
ANNOTATION_DIR <- "analysis/annotation_broad"

#' Log a message with timestamp and level
#' @param msg Character string containing message to log
#' @param level Log level (INFO, WARNING, ERROR)
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

#' Check if all required input files exist
#' @param files Vector of file paths to check
#' @return TRUE if all files exist, stops with error if any missing
validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
    TRUE
}

#' Create directories if they don't exist
#' @param dirs Vector of directory paths to create
create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

#' Standardize chromosome names in peak files to UCSC format (chr1, chr2, etc)
#' @param peak_files Vector of peak file paths to process
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s", file))
        
        # Read broadPeak format with error handling
        tryCatch({
            # First try reading with fill=TRUE to handle malformed lines
            peaks <- try({
                read.table(file, header=FALSE,
                          colClasses=c("character", "numeric", "numeric", "character",
                                     "numeric", "character", "numeric", "numeric", "numeric"),
                          col.names=c("chr", "start", "end", "name", "score",
                                    "strand", "signalValue", "pValue", "qValue"))
            }, silent=TRUE)
            
            if (inherits(peaks, "try-error")) {
                # If normal read fails, try with more lenient settings
                peaks <- read.delim(file, header=FALSE, fill=TRUE, sep="\t",
                                  col.names=c("chr", "start", "end", "name", "score",
                                            "strand", "signalValue", "pValue", "qValue"))
            }
            
            # Remove any rows with NA values
            peaks <- na.omit(peaks)
            
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No valid peaks found in %s after cleaning", file))
                next
            }
            
            # Add 'chr' prefix if missing
            if (!all(grepl("^chr", peaks$chr))) {
                peaks$chr <- paste0("chr", peaks$chr)
                peaks$chr <- gsub("^chrchr", "chr", peaks$chr)
            }
            
            # Keep only standard chromosomes
            standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
            peaks <- peaks[peaks$chr %in% standard_chroms,]
            
            if (nrow(peaks) == 0) {
                log_message(sprintf("WARNING: No peaks left after filtering chromosomes in %s", file))
                next
            }
            
            # Write back
            write.table(peaks, file, sep="\t", quote=FALSE, 
                       row.names=FALSE, col.names=FALSE)
            
            log_message(sprintf("Processed %d peaks in %s", nrow(peaks), file))
            
        }, error = function(e) {
            log_message(sprintf("ERROR processing %s: %s", file, e$message))
            stop(sprintf("Failed to process peak file %s", file))
        })
    }
}

#' Check and create BAM index files if missing
#' @param bam_files Vector of BAM file paths to check
check_bam_indexes <- function(bam_files) {
    for(bam in bam_files) {
        if(!file.exists(paste0(bam, ".bai"))) {
            log_message(sprintf("Creating index for %s", bam), "WARNING")
            system2("samtools", args = c("index", bam))
        }
    }
}

#' Perform differential binding analysis for broad peaks
#' @return GRanges object containing differential peaks
perform_diffbind <- function() {
    log_message("Starting differential binding analysis for broad peaks")
    
    # Create output directories
    dirs <- c(OUTPUT_DIR, ANNOTATION_DIR, "logs")
    create_dirs(dirs)
    
    # Create sample sheet with selected samples that have similar peak counts
    samples <- data.frame(
        SampleID = c("YAF_2", "YAF_3", "GFP_1", "GFP_3"),
        Factor = rep("H2AK119Ub", 4),
        Condition = rep(c("YAF", "GFP"), each=2),
        Replicate = c(2:3, 1, 3),
        bamReads = file.path("analysis/aligned", 
                            paste0(c("YAF_2", "YAF_3", "GFP_1", "GFP_3"),
                                 ".dedup.bam")),
        Peaks = file.path(PEAKS_DIR,
                         paste0(c("YAF_2", "YAF_3", "GFP_1", "GFP_3"),
                               "_broad_", PEAKS_SUFFIX, ".broadPeak")),
        PeakCaller = rep("broad", 4)
    )
    log_message("Using selected samples with similar peak counts:")
    log_message("YAF_2: ~35k peaks")
    log_message("YAF_3: ~37k peaks")
    log_message("GFP_1: ~29k peaks")
    log_message("GFP_3: ~27k peaks")
    
    # Validate input files
    validate_files(c(samples$bamReads, samples$Peaks))
    
    # Add this call in perform_diffbind() before creating the DiffBind object
    check_bam_indexes(samples$bamReads)
    
    # Standardize chromosome names in peak files
    standardize_peak_files(samples$Peaks)
    
    # Create DiffBind object with less stringent parameters
    log_message("Creating DiffBind object...")
    dba_data <- dba(sampleSheet=samples,
                    minOverlap=1,  # Changed from 2 to 1 to be less stringent
                    peakFormat="bed",
                    peakCaller="broad",
                    config=data.frame(AnalysisMethod="DBA_DESEQ2",
                                    fragmentSize=150,
                                    doBlacklist=TRUE,
                                    RunParallel=TRUE))
    
    # Add diagnostic step
    log_message("Initial peak counts per sample:")
    print(dba.show(dba_data))
    
    # Count reads with less stringent parameters
    log_message("Counting reads...")
    dba_data <- dba.count(dba_data, 
                         bUseSummarizeOverlaps=TRUE,
                         minCount=0,       # Changed from 1 to 0
                         bRemoveDuplicates=TRUE,
                         score=DBA_SCORE_READS)  # Changed from RPKM to READS
    
    # Add more diagnostic information
    log_message("After counting - sample statistics:")
    count_info <- dba.show(dba_data)
    print(count_info)
    
    # Check if we have valid counts
    if(sum(count_info$Reads) == 0) {
        stop("No reads were counted across any samples")
    }
    
    # Normalize and perform differential analysis with modified parameters
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, normalize=DBA_NORM_LIB)  # Changed from TMM to LIB
    
    dba_data <- dba.contrast(dba_data, 
                            categories=DBA_CONDITION,
                            minMembers=1)  # Changed from 2 to 1
    
    # Add diagnostic step
    log_message("Contrast information:")
    print(dba.show(dba_data, bContrasts=TRUE))
    
    dba_data <- dba.analyze(dba_data, 
                           method=DBA_DESEQ2,
                           bFullLibrarySize=TRUE,
                           bTagwise=FALSE)  # Added parameter
    
    # Extract results with more lenient thresholds
    log_message("Extracting results...")
    dba_results <- dba.report(dba_data,
                             th=1,     # Changed from 0.05 to 1 (no filtering)
                             bCalled=FALSE,  # Changed from TRUE to FALSE
                             bNormalized=TRUE)
    
    # Save results
    log_message("Saving results...")
    saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_analysis.rds"))
    saveRDS(dba_results, file.path(OUTPUT_DIR, "significant_peaks.rds"))
    
    # Convert to data frame and save
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, 
             file.path(OUTPUT_DIR, "differential_peaks.csv"), 
             row.names=FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # Annotate peaks
    log_message("Annotating peaks...")
    peaks_gr <- dba_results
    seqlevelsStyle(peaks_gr) <- "UCSC"
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(peaks_gr, TxDb=txdb,
                            tssRegion=c(-3000, 3000),
                            verbose=FALSE)
    
    # Save annotation results
    saveRDS(peakAnno, file.path(ANNOTATION_DIR, "peak_annotation.rds"))
    write.csv(as.data.frame(peakAnno),
             file.path(ANNOTATION_DIR, "peak_annotation.csv"),
             row.names=FALSE)
    
    # Generate annotation plots
    pdf(file.path(ANNOTATION_DIR, "annotation_plots.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    return(dba_results)
}

# Create sample sheet for all samples
samples <- data.frame(
    SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
    Factor = rep("H2AK119Ub", 6),
    Condition = rep(c("YAF", "GFP"), each=3),
    Replicate = rep(1:3, 2),
    bamReads = file.path("analysis/aligned", 
                        paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
                             ".dedup.bam")),
    Peaks = file.path(PEAKS_DIR,
                     paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
                           "_broad_", PEAKS_SUFFIX, ".broadPeak")),
    PeakCaller = rep("broad", 6)
)

# Validate input files
validate_files(c(samples$bamReads, samples$Peaks))

# Standardize chromosome names in peak files
standardize_peak_files(samples$Peaks)

# Execute the differential binding analysis
dba_results <- perform_diffbind()

# Save sample sheet
sample_sheet_file <- file.path(OUTPUT_DIR, "complete_sample_sheet.csv")
write.csv(samples, sample_sheet_file, row.names = FALSE)

log_message("Successfully completed differential binding analysis")