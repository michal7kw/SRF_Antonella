#' This script performs differential binding analysis on H2AK119Ub ChIP-seq data
#' comparing YAF and GFP samples. It identifies regions with significantly different
#' H2AK119Ub levels between conditions.
#'
#' Input files:
#' - analysis/5_peak_calling_strict/{sample}_broad_peaks_final.broadPeak: Final peak calls for each sample
#' - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample
#'
#' Output files:
#' - analysis/7_differential_binding_strict/: Directory containing DiffBind results
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
    PeakCaller = rep("macs", length(sample_ids))
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

# Thoroughly clean and fix peak files
fix_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("DEBUG: fix_peak_files: Processing file: %s", file))
        log_message(sprintf("DEBUG: fix_peak_files: File %s exists: %s", file, file.exists(file)))
        if (file.exists(file)) {
            log_message(sprintf("DEBUG: fix_peak_files: File %s is empty: %s", file, file.info(file)$size == 0))
        }
        log_message(sprintf("Fixing and standardizing peak file: %s", file))
        
        # Read the file line by line
        con <- file(file, "r")
        lines <- readLines(con)
        close(con)
        
        # Process each line to ensure proper formatting and type conversion
        # Initialize lists to store columns with correct types
        col_chr <- character(0)
        col_start <- integer(0)
        col_end <- integer(0)
        col_name <- character(0)
        col_score <- integer(0)
        col_strand <- character(0)
        col_signalValue <- numeric(0)
        col_pValue <- numeric(0)
        col_qValue <- numeric(0)

        for (line in lines) {
            # Skip empty lines
            if (nchar(trimws(line)) == 0) {
                next
            }
            
            fields <- strsplit(line, "\\s+")[[1]]
            
            if (length(fields) != 9) {
                if (length(fields) > 9 && grepl("^chr\\d+\\d+", fields[1])) {
                    chr_val <- sub("^(chr\\d+).*", "\\1", fields[1])
                    rest_val <- sub("^chr\\d+(\\d+).*", "\\1", fields[1])
                    fields <- c(chr_val, rest_val, fields[2:length(fields)])
                }
                if (length(fields) != 9) {
                    log_message(sprintf("Skipping malformed line in %s: %s", file, line), "WARNING")
                    next
                }
            }
            
            current_chr <- fields[1]
            if (grepl("^chr", current_chr)) {
                current_chr <- sub("^chr", "", current_chr)
            }
            if (is.na(current_chr) || nchar(trimws(current_chr)) == 0) {
                log_message(sprintf("Skipping line in %s due to empty or NA chromosome name (original line: '%s')", file, line), "WARNING")
                next
            }
            valid_chrom_pattern <- "^([1-9]|1[0-9]|2[0-2]|[XY]|M|MT)$"
            if (!grepl(valid_chrom_pattern, current_chr)) {
                log_message(sprintf("Skipping line in %s due to non-standard chromosome name '%s' (original line: '%s')", file, current_chr, line), "WARNING")
                next
            }
            col_chr <- c(col_chr, current_chr)
            col_name <- c(col_name, fields[4])
            col_strand <- c(col_strand, ifelse(fields[6] %in% c(".", "+", "-"), fields[6], "."))

            # Integer columns (start, end, score)
            start_val <- suppressWarnings(as.integer(fields[2])); if(is.na(start_val) || !is.finite(start_val)) start_val <- 0
            col_start <- c(col_start, start_val)
            end_val <- suppressWarnings(as.integer(fields[3])); if(is.na(end_val) || !is.finite(end_val)) end_val <- 0
            col_end <- c(col_end, end_val)
            score_val <- suppressWarnings(as.integer(fields[5])); if(is.na(score_val) || !is.finite(score_val)) score_val <- 0
            col_score <- c(col_score, score_val)

            # Numeric columns (signalValue, pValue, qValue)
            signal_val <- suppressWarnings(as.numeric(fields[7])); if(is.na(signal_val) || !is.finite(signal_val)) signal_val <- 0.0
            col_signalValue <- c(col_signalValue, signal_val)
            pval_val <- suppressWarnings(as.numeric(fields[8])); if(is.na(pval_val) || !is.finite(pval_val)) pval_val <- 0.0
            col_pValue <- c(col_pValue, pval_val)
            qval_val <- suppressWarnings(as.numeric(fields[9])); if(is.na(qval_val) || !is.finite(qval_val)) qval_val <- 0.0
            col_qValue <- c(col_qValue, qval_val)
        }
        
        # Write fixed data back to file using write.table
        log_message(sprintf("DEBUG: fix_peak_files: Before writing to %s, number of processed valid lines: %d", file, length(col_chr)))
        if (length(col_chr) > 0) {
            fixed_df <- data.frame(
                chr = col_chr,
                start = col_start,
                end = col_end,
                name = col_name,
                score = col_score,
                strand = col_strand,
                signalValue = col_signalValue,
                pValue = col_pValue,
                qValue = col_qValue,
                stringsAsFactors = FALSE
            )
            log_message(sprintf("DEBUG: fix_peak_files: First 3 rows of data.frame for %s:", file))
            print(head(fixed_df, 3))
            
            # Use write.table to ensure R's standard numeric->string conversion for file output
            write.table(fixed_df, file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
            
            log_message(sprintf("DEBUG: fix_peak_files: Successfully wrote data.frame to %s using write.table", file))
            log_message(sprintf("Successfully fixed %d lines in %s", nrow(fixed_df), file))
        } else {
            log_message(sprintf("DEBUG: fix_peak_files: No valid lines to write for %s", file))
            log_message(sprintf("ERROR: No valid lines found in %s after fixing", file), "ERROR")
            stop(sprintf("No valid peaks in %s", file))
        }
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

# Perform differential binding analysis for broad peaks with less stringent parameters
perform_diffbind <- function(samples) {
    log_message("Starting differential binding analysis for broad peaks with adjusted parameters")
    
    # Create output directories
    # dirs <- c(OUTPUT_DIR, ANNOTATION_DIR)
    # create_dirs(dirs)
    
    # Validate input files
    validate_files(c(samples$bamReads, samples$Peaks))
    
    # Check BAM indexes
    check_bam_indexes(samples$bamReads)
    
    # Fix peak files - more thorough cleaning for strict analysis
    fix_peak_files(samples$Peaks)
    
    # Create DiffBind object with less stringent parameters for strict analysis
    log_message("Creating DiffBind object with adjusted parameters...")
    log_message("DEBUG: perform_diffbind: Samples data frame before dba() call:")
    print(samples)
    log_message("DEBUG: perform_diffbind: Checking peak files listed in samples$Peaks before dba() call (content on disk after fix_peak_files):")
    for (peak_file_path_idx in 1:nrow(samples)) {
        peak_file_path <- samples$Peaks[peak_file_path_idx]
        log_message(sprintf("DEBUG: perform_diffbind: Peak file path: %s", peak_file_path))
        log_message(sprintf("DEBUG: perform_diffbind: Exists: %s", file.exists(peak_file_path)))
        if (file.exists(peak_file_path)) {
            log_message(sprintf("DEBUG: perform_diffbind: Is empty: %s", file.info(peak_file_path)$size == 0))
            if (file.info(peak_file_path)$size > 0) {
                log_message(sprintf("DEBUG: perform_diffbind: First 5 lines of %s on disk:", peak_file_path))
                print(readLines(peak_file_path, n=5))
            }
        }
    }
    
    dba_data <- dba(sampleSheet = samples,
                    minOverlap = 1,    # Keep this lenient to include more sites
                    peakFormat = "broadPeak",
                    peakCaller = "macs",
                    scoreCol = 7, # Explicitly use signalValue (column 7) for broadPeak
                    config = data.frame(AnalysisMethod = "DBA_DESEQ2",
                                        bUseStandardChromosomes = FALSE, # Tell DiffBind to use chrom names as is
                                        fragmentSize = 150,
                                        doBlacklist = TRUE,
                                        RunParallel = TRUE,
                                        summits = FALSE))  # Disable summit calculation for broad peaks
    
    log_message("DEBUG: perform_diffbind: dba_data object created.")
    log_message("DEBUG: perform_diffbind: dba_data$config:")
    print(dba_data$config)
    
    if (!is.null(dba_data$peaks)) {
        log_message("DEBUG: perform_diffbind: Inspecting dba_data$peaks:")
        for (i in seq_along(dba_data$peaks)) {
            log_message(sprintf("DEBUG: perform_diffbind: Sample %d: %s", i, names(dba_data$peaks)[i]))
            if (!is.null(dba_data$peaks[[i]]$peaks)) {
                log_message(sprintf("DEBUG: perform_diffbind: Names of peak data for sample %d: %s", i, paste(names(dba_data$peaks[[i]]), collapse=", ")))
                # Attempt to get seqlevels, handle potential NULL or non-GRanges object
                current_peaks <- dba_data$peaks[[i]]$peaks
                if (inherits(current_peaks, "GRanges")) {
                     log_message(sprintf("DEBUG: perform_diffbind: Seqlevels for sample %d: %s", i, paste(seqlevels(current_peaks), collapse=", ")))
                } else if (!is.null(current_peaks)) {
                    log_message(sprintf("DEBUG: perform_diffbind: dba_data$peaks[[%d]]$peaks is not a GRanges object, it is a %s. Cannot get seqlevels directly.", i, class(current_peaks)))
                } else {
                    log_message(sprintf("DEBUG: perform_diffbind: dba_data$peaks[[%d]]$peaks is NULL.", i))
                }
            } else {
                log_message(sprintf("DEBUG: perform_diffbind: dba_data$peaks[[%d]]$peaks is NULL.", i))
            }
        }
    } else {
        log_message("DEBUG: perform_diffbind: dba_data$peaks is NULL.")
    }
    log_message("DEBUG: perform_diffbind: Attempting direct import of peak files using rtracklayer::import():", "DEBUG")
    for (i in 1:nrow(samples)) {
        peak_file_path <- samples$Peaks[i]
        log_message(sprintf("DEBUG: perform_diffbind: Processing peak file with rtracklayer: %s", peak_file_path), "DEBUG")
        tryCatch({
            # Ensure the file exists and is not empty before attempting import
            if (!file.exists(peak_file_path)) {
                log_message(paste("DEBUG: ERROR: Peak file does not exist for rtracklayer::import():", peak_file_path), "ERROR")
                next # Skip to the next iteration
            }
            if (file.info(peak_file_path)$size == 0) {
                log_message(paste("DEBUG: ERROR: Peak file is empty for rtracklayer::import():", peak_file_path), "ERROR")
                next # Skip to the next iteration
            }
            
            # Use format = "broadPeak" for diagnostic import, reflecting dba() parameters
            gr <- rtracklayer::import(peak_file_path, format = "broadPeak")
            log_message(paste("DEBUG: Successfully imported with rtracklayer::import(format='broadPeak'):", peak_file_path), "DEBUG")
            log_message("DEBUG: seqinfo from rtracklayer::import(format='broadPeak'):", "DEBUG")
            print(seqinfo(gr))
            log_message("DEBUG: head(gr, 3) from rtracklayer::import(format='broadPeak'):", "DEBUG")
            print(head(gr, 3))
        }, error = function(e) {
            log_message(paste("DEBUG: ERROR importing file with rtracklayer::import(format='broadPeak'):", peak_file_path, "-", e$message), "DEBUG")
        })
    }
    
    log_message("Initial peak counts per sample:")
    print(dba.show(dba_data))
    
    # Count reads with less stringent parameters
    log_message("Counting reads with adjusted parameters...")
    dba_data <- dba.count(dba_data,
                          bUseSummarizeOverlaps = TRUE,
                          filter = 0, # Explicitly set filter
                          minCount = 0,       # Keep for clarity
                          bRemoveDuplicates = TRUE,
                          score = DBA_SCORE_READS,
                          summits = FALSE)
   
   log_message("DEBUG: perform_diffbind: dba.count finished.")
   log_message("DEBUG: perform_diffbind: dba_data object summary after dba.count:")
   print(dba_data)
   
   if (!is.null(dba_data$merged) && inherits(dba_data$merged, "GRanges")) {
       log_message("DEBUG: perform_diffbind: Seqlevels of dba_data$merged after dba.count:")
       print(seqlevels(dba_data$merged))
   } else {
       log_message("DEBUG: perform_diffbind: dba_data$merged is NULL or not a GRanges object after dba.count.")
   }
   
   if (!is.null(dba_data$binding)) {
       log_message("DEBUG: perform_diffbind: Dimensions of dba_data$binding after dba.count:")
       print(dim(dba_data$binding))
       log_message("DEBUG: perform_diffbind: Summary of dba_data$binding after dba.count:")
       print(summary(dba_data$binding))
   } else {
       log_message("DEBUG: perform_diffbind: dba_data$binding is NULL after dba.count.")
   }
   
   log_message("After counting - sample statistics:")
    count_info <- dba.show(dba_data)
    print(count_info)
    
    if(sum(count_info$Reads) == 0) {
        stop("No reads were counted across any samples")
    }
    log_message("Full dba_data object state after dba.count:")
    print(dba_data)
    
    # Normalize and perform differential analysis with less stringent parameters
    log_message("Performing differential analysis with adjusted parameters...")
    dba_data <- dba.normalize(dba_data, 
                              normalize = DBA_NORM_LIB,
                              background = FALSE)  # Disable background normalization for strict analysis
    
    log_message("Setting up contrasts...")
    dba_data <- dba.contrast(dba_data, 
                             categories = DBA_CONDITION,
                             minMembers = 2)
    
    log_message("Contrast information:")
    print(dba.show(dba_data, bContrasts = TRUE))

    # Add this before dba.analyze
    log_message("Binding site statistics before dba.analyze (summary of DBA object):")
    print(summary(dba_data)) # Changed from dba.count(dba_data)

    # Save a baseline count report for debugging
    write.csv(dba.peakset(dba_data, bRetrieve=TRUE), 
            file.path(OUTPUT_DIR, "pre_analysis_binding_sites.csv"))
    
    log_message("Inspecting dba_data object just before primary dba.analyze call:")
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
    
    tryCatch({
        log_message("Starting DESeq2 analysis with default DiffBind settings...")
        
        # Before analysis, inspect count distribution from the DBA object
        # This reflects the matrix that will be used by dba.analyze
        if(!is.null(dba_data$binding) && nrow(dba_data$binding) > 0 && ncol(dba_data$binding) > 0) {
            log_message(sprintf("Internal count matrix (dba_data$binding) dimensions: %d sites x %d samples",
                                nrow(dba_data$binding), ncol(dba_data$binding)))
            # Save this internal count matrix for troubleshooting if needed
            write.csv(as.data.frame(dba_data$binding), file.path(OUTPUT_DIR, "internal_raw_count_matrix_before_analyze.csv"), row.names = TRUE)
        } else {
            log_message("WARNING: dba_data$binding is NULL or empty before dba.analyze call.")
        }
        
        log_message("DEBUG: perform_diffbind: Preparing for dba.analyze.")
        if (!is.null(dba_data$config$doBlacklist) && dba_data$config$doBlacklist) {
            log_message("DEBUG: perform_diffbind: Blacklisting is enabled (dba_data$config$doBlacklist == TRUE).")
            if (!is.null(dba_data$blacklist) && !is.null(dba_data$blacklist$regions) && inherits(dba_data$blacklist$regions, "GRanges")) {
                log_message("DEBUG: perform_diffbind: Seqlevels of dba_data$blacklist$regions:")
                print(seqlevels(dba_data$blacklist$regions))
            } else {
                log_message("DEBUG: perform_diffbind: dba_data$blacklist$regions is NULL or not a GRanges object.")
            }
        } else {
            log_message("DEBUG: perform_diffbind: Blacklisting is NOT enabled or dba_data$config$doBlacklist is NULL.")
        }
        
        if (!is.null(dba_data$merged) && inherits(dba_data$merged, "GRanges")) {
            log_message("DEBUG: perform_diffbind: Seqlevels of dba_data$merged before dba.analyze:")
            print(seqlevels(dba_data$merged))
        } else {
            log_message("DEBUG: perform_diffbind: dba_data$merged is NULL or not a GRanges object before dba.analyze.")
        }

        # Simplest dba.analyze call for DESeq2
        # This allows DESeq2 to use its own default filtering mechanisms
        dba_data_analyzed <- dba.analyze(dba_data, method = DBA_DESEQ2)
        
        log_message("DESeq2 analysis (dba.analyze) completed successfully.")
        
        # Get detailed statistics about the analysis
        analysis_stats <- dba.show(dba_data_analyzed, bContrasts = TRUE)
        log_message("Analysis statistics:")
        print(analysis_stats)
        
        # Plot dispersion estimates using the analyzed object
        pdf(file.path(OUTPUT_DIR, "differential_analysis_plots.pdf"))
        tryCatch({
            log_message("Plotting Heatmap...")
            dba.plotHeatmap(dba_data_analyzed, contrast=1) # Plot for the first contrast
            title("Sample Correlation Heatmap (Post-Analysis)")
            
            log_message("Plotting PCA...")
            dba.plotPCA(dba_data_analyzed, contrast=1, label=DBA_ID)
            title("PCA Plot (Post-Analysis)")
            
            log_message("Plotting MA Plot...")
            dba.plotMA(dba_data_analyzed, contrast=1)
            title("MA Plot (Post-Analysis)")
            
            log_message("Plotting Volcano Plot...")
            dba.plotVolcano(dba_data_analyzed, contrast=1)
            title("Volcano Plot (Post-Analysis)")
            
        }, error = function(e_plot) {
            log_message(sprintf("Warning: Could not generate some post-analysis plots: %s", e_plot$message), "WARNING")
        }, finally = {
            if(length(dev.list()) > 0) dev.off() # Ensure device is closed if plots were attempted
        })
        
        # Assign the analyzed object back to dba_data for reporting
        dba_data <- dba_data_analyzed
        
    }, error = function(e_analyze) {
        log_message(sprintf("Error in differential analysis (dba.analyze): %s", e_analyze$message), "ERROR")
        # Save the DBA object *before* the error for inspection
        saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_object_before_analyze_error.rds"))
        log_message("Saved dba_data object (pre-analysis attempt) to dba_object_before_analyze_error.rds for inspection.")
        stop(sprintf("Differential analysis (dba.analyze) failed: %s", e_analyze$message))
    })
    
    log_message("Extracting results with less stringent threshold...")
    # Use a more lenient threshold for the strict analysis
    dba_results <- dba.report(dba_data,
                             th = 1,  # No p-value threshold
                             bCalled = FALSE,
                             bNormalized = TRUE,
                             bCounts = TRUE)  # Include read counts
    
    log_message("Saving results...")
    saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_analysis.rds"))
    saveRDS(dba_results, file.path(OUTPUT_DIR, "significant_peaks.rds"))
    
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, file.path(OUTPUT_DIR, "differential_peaks.csv"), row.names = FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # === Peak Annotation ===
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
