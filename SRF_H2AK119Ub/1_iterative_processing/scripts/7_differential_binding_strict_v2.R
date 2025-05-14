
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

# Thoroughly clean and fix peak files
fix_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Fixing and standardizing peak file: %s", file))
        
        # Read the file line by line
        con <- file(file, "r")
        lines <- readLines(con)
        close(con)
        
        # Process each line to ensure proper formatting
        fixed_lines <- character(0)
        for (line in lines) {
            # Skip empty lines
            if (nchar(trimws(line)) == 0) {
                next
            }
            
            # Split by any whitespace (handles both spaces and tabs)
            fields <- strsplit(line, "\\s+")[[1]]
            
            # Ensure we have exactly 9 fields for broadPeak format
            if (length(fields) != 9) {
                # If we have more than 9 fields, try to fix common issues
                if (length(fields) > 9) {
                    # Check if chromosome field might be malformed (e.g., "chr13754338089")
                    if (grepl("^chr\\d+\\d+", fields[1])) {
                        # Extract chromosome number
                        chr <- sub("^(chr\\d+).*", "\\1", fields[1])
                        # Extract the rest as potential start position
                        rest <- sub("^chr\\d+(\\d+).*", "\\1", fields[1])
                        
                        # Reconstruct fields with proper separation
                        fields <- c(chr, rest, fields[2:length(fields)])
                    }
                }
                
                # If we still don't have 9 fields, skip this line
                if (length(fields) != 9) {
                    log_message(sprintf("Skipping malformed line in %s: %s", file, line), "WARNING")
                    next
                }
            }
            
            # Standardize chromosome format: REMOVE "chr" prefix
            # This assumes BAM files use "1", "2", "MT" etc.
            if (grepl("^chr", fields[1])) {
                fields[1] <- sub("^chr", "", fields[1])
            }
            
            # Ensure chromosome name is not empty or NA after 'chr' removal or originally
            if (is.na(fields[1]) || nchar(trimws(fields[1])) == 0) {
                log_message(sprintf("Skipping line in %s due to empty or NA chromosome name (original line: '%s')", file, line), "WARNING")
                next
            }

            # Add stricter validation for chromosome names (e.g., for hg38 primary assembly)
            # This regex matches 1-22, X, Y, M (or MT).
            valid_chrom_pattern <- "^([1-9]|1[0-9]|2[0-2]|[XY]|M|MT)$"
            if (!grepl(valid_chrom_pattern, fields[1])) {
                log_message(sprintf("Skipping line in %s due to non-standard chromosome name '%s' after 'chr' removal (original line: '%s')", file, fields[1], line), "WARNING")
                next
            }
            
            # Ensure numeric fields are valid numbers
            for (i in c(2, 3, 5, 7, 8, 9)) {
                if (!grepl("^[0-9.]+$", fields[i])) {
                    fields[i] <- "0"  # Replace invalid numbers with 0
                }
            }
            
            # Ensure strand is valid
            if (!fields[6] %in% c(".", "+", "-")) {
                fields[6] <- "."
            }
            
            # Join fields with tabs
            fixed_line <- paste(fields, collapse="\t")
            fixed_lines <- c(fixed_lines, fixed_line)
        }
        
        # Write fixed lines back to file
        if (length(fixed_lines) > 0) {
            writeLines(fixed_lines, file)
            log_message(sprintf("Successfully fixed %d lines in %s", length(fixed_lines), file))
        } else {
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
    dba_data <- dba(sampleSheet = samples,
                    minOverlap = 1,    # Keep this lenient to include more sites
                    peakFormat = "bed",      # 'bed' is fine as MACS broadPeak is BED-like
                    peakCaller = "macs",       # Corrected to 'macs' for broadPeak files
                    config = data.frame(AnalysisMethod = "DBA_DESEQ2",
                                        bUseStandardChromosomes = FALSE, # Tell DiffBind to use chrom names as is
                                        fragmentSize = 150,
                                        doBlacklist = TRUE,
                                        RunParallel = TRUE,
                                        summits = FALSE))  # Disable summit calculation for broad peaks
    
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

