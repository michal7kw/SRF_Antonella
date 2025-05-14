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

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("Usage: Rscript script.R OUTPUT_DIR PEAKS_DIR PEAKS_SUFFIX ALIGN_DIR SAMPLE_IDS SAMPLE_CONDITIONS SAMPLE_REPLICATES")
}

# Define global constants from command line arguments
OUTPUT_DIR <- args[1]
PEAKS_DIR <- args[2]
PEAKS_SUFFIX <- args[3]
ALIGN_DIR <- args[4]
# ANNOTATION_DIR <- file.path(OUTPUT_DIR, "./annotation")

# Parse sample information
sample_ids <- strsplit(args[5], ",")[[1]]
sample_conditions <- strsplit(args[6], ",")[[1]]
sample_replicates <- as.numeric(strsplit(args[7], ",")[[1]])

if (length(sample_ids) != length(sample_conditions) || length(sample_ids) != length(sample_replicates)) {
    stop("Sample IDs, conditions, and replicates must have the same length")
}

# Create samples data frame
samples <- data.frame(
    SampleID = sample_ids,
    Factor = rep("H2AK119Ub", length(sample_ids)),
    Condition = sample_conditions,
    Replicate = sample_replicates,
    bamReads = file.path(ALIGN_DIR, paste0(sample_ids, ".dedup.bam")),
    Peaks = file.path(PEAKS_DIR, paste0(sample_ids, "_broad_", PEAKS_SUFFIX, ".broadPeak")),
    PeakCaller = rep("broad", length(sample_ids))
)

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
# This function will standardize chromosome names to non-UCSC style (e.g., "1", "X", "MT").
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s to non-UCSC style (e.g., 1, X, MT)", file))
        tryCatch({
            peaks <- NULL
            error_lines <- c()
            
            # Read the file line by line first to identify problematic lines
            con <- file(file, "r")
            lines <- readLines(con)
            close(con)
            
            valid_lines <- character(0)
            for (i in seq_along(lines)) {
                line <- lines[i]
                fields <- strsplit(line, "\t")[[1]]
                
                if (length(fields) != 9) {
                    error_lines <- c(error_lines,
                        sprintf("Line %d: Expected 9 fields, found %d", i, length(fields)))
                    next
                }
                
                # Validate chromosome format (allow 'chr1' or '1' as input)
                # We will strip "chr" for validation and later for standardization.
                chr_value_for_validation <- gsub("^chr", "", fields[1])
                if (!grepl("^([1-9][0-9]?|X|Y|M|MT)$", chr_value_for_validation)) {
                    error_lines <- c(error_lines,
                        sprintf("Line %d: Invalid chromosome format: %s (validated as %s)", i, fields[1], chr_value_for_validation))
                    next
                }
                
                # Validate numeric fields
                if (!all(grepl("^[0-9]+$", fields[c(2,3,5)]) &
                        grepl("^[0-9.]+$", fields[c(7,8,9)]))) {
                    error_lines <- c(error_lines,
                        sprintf("Line %d: Invalid numeric fields", i))
                    next
                }
                
                valid_lines <- c(valid_lines, line)
            }
            
            if (length(error_lines) > 0) {
                error_file <- paste0(file, ".errors")
                writeLines(c(
                    sprintf("Errors found in %s:", file),
                    error_lines,
                    sprintf("\nTotal lines: %d", length(lines)),
                    sprintf("Valid lines: %d", length(valid_lines)),
                    sprintf("Invalid lines: %d", length(error_lines))
                ), error_file)
                
                log_message(sprintf("Found %d invalid lines in %s. See %s for details",
                    length(error_lines), file, error_file))
            }
            
            if (length(valid_lines) == 0) {
                stop(sprintf("No valid peaks found in %s after validation", file))
            }
            
            # Write valid lines to temporary file
            temp_file <- tempfile()
            writeLines(valid_lines, temp_file)
            
            # Read the cleaned data
            peaks_df <- read.table(temp_file, header=FALSE,
                colClasses=c("character", "numeric", "numeric", "character",
                            "numeric", "character", "numeric", "numeric", "numeric"),
                col.names=c("chr", "start", "end", "name", "score",
                            "strand", "signalValue", "pValue", "qValue"))
            
            unlink(temp_file)  # Clean up temporary file
            
            # Standardize chromosome names to non-UCSC style (e.g., "1", "X", "MT")
            original_nrow <- nrow(peaks_df)
            peaks_df$chr <- gsub("^chr", "", peaks_df$chr)
            log_message(sprintf("Removed 'chr' prefix for standardization in %s. Rows before: %d", file, original_nrow))

            # Filter for standard chromosomes (now in non-UCSC format)
            standard_chroms_non_ucsc <- c(as.character(1:22), "X", "Y", "MT")
            peaks_df <- peaks_df[peaks_df$chr %in% standard_chroms_non_ucsc, ]
            log_message(sprintf("Filtered from %d to %d peaks based on standard non-UCSC chromosomes in %s", original_nrow, nrow(peaks_df), file))
            
            if (nrow(peaks_df) == 0) {
                stop(sprintf("No peaks left after standardizing and filtering chromosomes in %s", file))
            }
            
            # Write back the cleaned data
            write.table(peaks_df, file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
            
            log_message(sprintf("Successfully processed and standardized %d peaks in %s to non-UCSC style", nrow(peaks_df), file))
            
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
    
    # # Create output directories
    # dirs <- c(OUTPUT_DIR, ANNOTATION_DIR)
    # create_dirs(dirs)
    
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
    tryCatch({
        dba_data <- dba.count(dba_data,
                              bUseSummarizeOverlaps = TRUE,
                              filter = 0, # Explicitly set filter
                              minCount = 0, # Keep for clarity, though filter=0 should cover it
                              bRemoveDuplicates = TRUE,
                              score = DBA_SCORE_READS,
                              summits = FALSE)
        
        log_message("dba.count completed. After counting - sample statistics:")
        count_info <- dba.show(dba_data)
        print(count_info)
        
        if(is.null(count_info) || sum(count_info$Reads) == 0) {
             log_message("WARNING: No reads were counted across any samples, or count_info is NULL.")
        }
        
        log_message("Full dba_data object state after dba.count:")
        print(dba_data)
        
        if (is.null(dba_data$binding) || nrow(dba_data$binding) == 0) {
            log_message("CRITICAL: dba_data$binding (counts table) is NULL or has 0 rows immediately after dba.count. This is the primary issue.")
            stop("CRITICAL: Counts table (dba_data$binding) is empty after dba.count.")
        }

    }, error = function(e) {
        log_message(sprintf("ERROR during dba.count: %s", e$message), level = "ERROR")
        log_message("Printing dba_data object state before dba.count error:", level = "ERROR")
        print(dba_data) # Print the state of dba_data before the error
        stop(sprintf("Failed during dba.count: %s", e$message))
    })
    
    # Normalize and perform differential analysis
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, 
                          normalize = DBA_NORM_LIB,
                          background = FALSE)  # Set to FALSE with no controls
    
    log_message("Setting up contrasts...")
    dba_data <- dba.contrast(dba_data, 
                             categories = DBA_CONDITION,
                             minMembers = 2)
    
    log_message("Contrast information:")
    print(dba.show(dba_data, bContrasts = TRUE))
    
    tryCatch({
        log_message("Inspecting dba_data object just before dba.analyze call:")
        if(is.null(dba_data$binding) || nrow(dba_data$binding) == 0) {
            log_message("ERROR: dba_data$binding (counts) is null or has 0 rows before dba.analyze!")
        } else {
            log_message(sprintf("Number of sites (rows in count table) in dba_data: %d", nrow(dba_data$binding)))
            log_message(sprintf("Number of samples (cols in count table) in dba_data: %d", ncol(dba_data$binding)))
        }
        if (!is.null(dba_data$DESeq2$object)) {
            log_message("DESeq2 object already exists in dba_data, which is unexpected here.")
        } else {
            log_message("DESeq2 object does not exist yet in dba_data (expected).")
        }
        
        log_message("Starting DESeq2 analysis with dispersion estimation...")
        
        # Analyze with DESeq2 method with correct parameters
        dba_data <- dba.analyze(dba_data,
                       method = DBA_DESEQ2)
        
        log_message("DESeq2 analysis completed successfully")
        
        # Get detailed statistics about the analysis
        analysis_stats <- dba.show(dba_data, bContrasts = TRUE)
        log_message("Analysis statistics:")
        print(analysis_stats)
        
        # Plot dispersion estimates
        pdf(file.path(OUTPUT_DIR, "dispersion_plots.pdf"))
        tryCatch({
            dba.plotHeatmap(dba_data)
            title("Sample correlation heatmap")
            
            dba.plotPCA(dba_data, DBA_CONDITION, label=DBA_ID)
            title("PCA plot of samples")
            
            dba.plotMA(dba_data)
            title("MA plot of differential binding")
            
            dba.plotVolcano(dba_data)
            title("Volcano plot of differential binding")
        }, error = function(e) {
            log_message(sprintf("Warning: Could not generate some plots: %s", e$message), "WARNING")
        }, finally = {
            dev.off()
        })
        
    }, error = function(e) {
        log_message(sprintf("Error in differential analysis: %s", e$message), "ERROR")
        stop(sprintf("Differential analysis failed: %s", e$message))
    })
    
    log_message("Extracting results...")
    dba_results <- dba.report(dba_data,
                         th = 1,
                         bCalled = FALSE,
                         bNormalized = TRUE,
                         bCounts = TRUE)  # Add this to include read counts
    
    log_message("Saving results...")
    saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_analysis.rds"))
    saveRDS(dba_results, file.path(OUTPUT_DIR, "significant_peaks.rds"))
    
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, file.path(OUTPUT_DIR, "differential_peaks.csv"), row.names = FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # # === Peak Annotation ===
    # log_message("Annotating peaks...")
    # peaks_gr <- dba_results
    
    # # For annotation we need UCSC-style names. Convert the non-UCSC (e.g. "1", "2", …)
    # # to UCSC style by adding the "chr" prefix.
    # seqlevels(peaks_gr) <- paste0("chr", seqlevels(peaks_gr))
    
    # txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    # peakAnno <- annotatePeak(peaks_gr, TxDb = txdb,
    #                          tssRegion = c(-3000, 3000),
    #                          verbose = FALSE)
    
    # saveRDS(peakAnno, file.path(ANNOTATION_DIR, "peak_annotation.rds"))
    # write.csv(as.data.frame(peakAnno),
    #           file.path(ANNOTATION_DIR, "peak_annotation.csv"),
    #           row.names = FALSE)
    
    # pdf(file.path(ANNOTATION_DIR, "annotation_plots.pdf"))
    # plotAnnoPie(peakAnno)
    # plotDistToTSS(peakAnno)
    # dev.off()
    
    return(dba_results)
}

# Execute the differential binding analysis
dba_results <- perform_diffbind(samples)

# Save sample sheet
sample_sheet_file <- file.path(OUTPUT_DIR, "complete_sample_sheet.csv")
write.csv(samples, sample_sheet_file, row.names = FALSE)

log_message("Successfully completed differential binding analysis")
