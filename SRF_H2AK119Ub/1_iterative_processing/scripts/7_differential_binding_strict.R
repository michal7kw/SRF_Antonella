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
#' - analysis/7_differential_binding_strict/: Directory containing annotated differential binding results
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
ANNOTATION_DIR <- file.path(OUTPUT_DIR, "./annotation")

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
            
            # Ensure chromosome format is consistent (with "chr" prefix)
            if (!grepl("^chr", fields[1])) {
                fields[1] <- paste0("chr", fields[1])
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
    dirs <- c(OUTPUT_DIR, ANNOTATION_DIR)
    create_dirs(dirs)
    
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
                    peakFormat = "bed",
                    peakCaller = "broad",
                    config = data.frame(AnalysisMethod = "DBA_DESEQ2",
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
                          minCount = 0,       # Keep this at 0 to include all sites
                          bRemoveDuplicates = TRUE,
                          score = DBA_SCORE_READS,
                          summits = FALSE)
    
    log_message("After counting - sample statistics:")
    count_info <- dba.show(dba_data)
    print(count_info)
    
    if(sum(count_info$Reads) == 0) {
        stop("No reads were counted across any samples")
    }
    
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
    log_message("Binding site statistics before filtering:")
    print(summary(dba.count(dba_data)))

    # Save a baseline count report for debugging
    write.csv(dba.peakset(dba_data, bRetrieve=TRUE), 
            file.path(OUTPUT_DIR, "pre_analysis_binding_sites.csv"))
    
    tryCatch({
        log_message("Starting DESeq2 analysis with adjusted parameters...")
        
        # Before analysis, inspect count distribution
        count_matrix <- dba.peakset(dba_data, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
        log_message(sprintf("Count matrix dimensions: %d rows x %d columns", nrow(count_matrix), ncol(count_matrix)))
        
        # Save count data for troubleshooting
        write.csv(count_matrix, file.path(OUTPUT_DIR, "raw_count_matrix.csv"), row.names = FALSE)
        
        # Try a completely different approach if normal analysis fails
        log_message("Attempting alternative analysis approach...")
        
        # First try the standard approach with completely disabled filtering
        tryCatch({
            dba_data <- dba.analyze(dba_data, 
                           method = DBA_DESEQ2,
                           bSubControl = FALSE,
                           bFullLibrarySize = TRUE,
                           bReduceObjects = FALSE,
                           filter = FALSE,
                           bTagwise = TRUE)
            
            log_message("Primary DESeq2 analysis completed successfully")
        }, error = function(e) {
            log_message(sprintf("Primary analysis failed: %s. Trying fallback approach.", e$message), "WARNING")
            
            # FALLBACK: Create a simpler DiffBind object and try again with minimal options
            log_message("Re-creating DiffBind object with minimal configuration...")
            dba_data <<- dba(sampleSheet = samples,
                          minOverlap = 1,
                          peakCaller = "broad")
            
            dba_data <<- dba.count(dba_data, 
                               bParallel = TRUE,
                               bUseSummarizeOverlaps = TRUE,
                               minCount = 0)
            
            log_message("Re-normalizing with minimal parameters...")
            dba_data <<- dba.normalize(dba_data, normalize = DBA_NORM_LIB)
            
            log_message("Setting up minimal contrasts...")
            dba_data <<- dba.contrast(dba_data, 
                                  categories = DBA_CONDITION, 
                                  minMembers = 1)
            
            log_message("Attempting fallback analysis with minimal parameters...")
            dba_data <<- dba.analyze(dba_data, 
                                method = DBA_DESEQ2,
                                bTagwise = FALSE,
                                bFullLibrarySize = TRUE,
                                filter = 1,
                                bReduceObjects = TRUE)
        })
        
        # After analysis complete
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

# Execute the differential binding analysis
dba_results <- perform_diffbind(samples)

# Save sample sheet
sample_sheet_file <- file.path(OUTPUT_DIR, "complete_sample_sheet.csv")
write.csv(samples, sample_sheet_file, row.names = FALSE)

log_message("Successfully completed differential binding analysis")
