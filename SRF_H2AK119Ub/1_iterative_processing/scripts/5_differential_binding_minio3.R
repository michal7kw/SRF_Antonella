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

# Constants
PEAKS_DIR <- "analysis/peaks2_improved"
PEAKS_SUFFIX <- "peaks_final"
OUTPUT_DIR <- "analysis/diffbind_broad"
ANNOTATION_DIR <- "analysis/annotation_broad"

# Utility functions
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
    TRUE
}

create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Function to standardize chromosome names in peak files
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s", file))
        
        # Read peak file with improved error handling and type specifications.
        peaks <- tryCatch({
            read.table(file, header = FALSE, stringsAsFactors = FALSE,
                       colClasses = c("character", "numeric", "numeric", "character",
                                      "numeric", "character", "numeric", "numeric", "numeric"),
                       col.names = c("chr", "start", "end", "name", "score",
                                     "strand", "foldChange", "pValue", "qValue"))
        }, error = function(e) {
            log_message(sprintf("Error reading file %s: %s", file, e$message), level = "ERROR")
            stop(e)
        })
        
        if (nrow(peaks) == 0) {
            log_message(sprintf("No peaks found in %s, skipping standardization.", file), level = "WARNING")
            next
        }
        
        # Add 'chr' prefix if missing and remove duplicate prefixes.
        peaks$chr <- ifelse(grepl("^chr", peaks$chr), peaks$chr, paste0("chr", peaks$chr))
        peaks$chr <- gsub("^chr+", "chr", peaks$chr)
        
        # Keep only standard chromosomes
        standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
        peaks <- peaks[peaks$chr %in% standard_chroms, ]
        
        # Reset row names
        rownames(peaks) <- NULL
        
        # Write back the standardized data to the same file
        write.table(peaks, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        log_message(sprintf("Processed %d peaks in %s", nrow(peaks), file))
    }
}

#' Perform differential binding analysis for broad peaks
perform_diffbind <- function() {
    log_message("Starting differential binding analysis for broad peaks")
    
    # Create output directories
    dirs <- c(OUTPUT_DIR, ANNOTATION_DIR, "logs")
    create_dirs(dirs)
    
    # Create sample sheet
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
    
    # Create DiffBind object with optimized parameters for broad peaks
    log_message("Creating DiffBind object...")
    dba_data <- dba(sampleSheet=samples,
                    minOverlap=2,  # Require peaks to be present in at least 2 samples
                    peakFormat="bed",
                    peakCaller="broad",
                    config=data.frame(AnalysisMethod="DBA_DESEQ2",
                                    fragmentSize=150,
                                    doBlacklist=TRUE,
                                    RunParallel=TRUE))
    
    # Count reads with parameters suitable for broad peaks
    log_message("Counting reads...")
    dba_data <- dba.count(dba_data, 
                         bUseSummarizeOverlaps=TRUE,
                         minCount=1,
                         bRemoveDuplicates=TRUE,
                         score=DBA_SCORE_RPKM)  # Use RPKM for broad peaks
    
    # Print diagnostic information
    count_info <- dba.show(dba_data)
    log_message("Read counts and FRiP per sample:")
    print(count_info[,c("Reads", "FRiP")])
    
    # Normalize and perform differential analysis
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, normalize=DBA_NORM_TMM)
    
    dba_data <- dba.contrast(dba_data, 
                            categories=DBA_CONDITION,
                            minMembers=2)
    
    dba_data <- dba.analyze(dba_data, 
                           method=DBA_DESEQ2,
                           bFullLibrarySize=TRUE)
    
    # Extract results with FDR < 0.05 and absolute fold change > 1.5
    log_message("Extracting significant results...")
    dba_results <- dba.report(dba_data,
                             th=0.05,     # FDR threshold
                             bCalled=TRUE,
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

# Run analysis
results <- perform_diffbind()