#!/usr/bin/env Rscript

#####################################################################
# Cross-reference Analysis Script for V5 and H2AK119Ub ChIP-seq Data
#####################################################################

# DESCRIPTION:
# This script performs cross-reference analysis between V5 ChIP-seq 
# and H2AK119Ub ChIP-seq data to identify overlapping binding sites

# INPUT FILES:
# - V5 ChIP-seq peaks:
#   * ../SRF_V5/results/peaks/SES-V5ChIP-Seq2_S6_narrow_peaks.narrowPeak
#   * ../SRF_V5/results/peaks/SES-V5ChIP-Seq2_S6_broad_peaks.broadPeak
#
# - V5 ChIP-seq signal:
#   * ../SRF_V5/results/bigwig/SES-V5ChIP-Seq2_S6.bw (ChIP)
#   * ../SRF_V5/results/bigwig/InputSES-V5ChIP-Seq_S2.bw (Input)
#
# - H2AK119Ub peaks (broad):
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/peaks/GFP_[1-3]_broad_peaks_final.broadPeak
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/peaks/YAF_[1-3]_broad_peaks_final.broadPeak
#
# - H2AK119Ub signal:
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/visualization/GFP_[1-3].bw
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/visualization/YAF_[1-3].bw

# OUTPUT FILES:
# - Processed data:
#   * results/processed_peaks.rds
#   * results/normalized_v5_signal.rds
#   * results/*_categorized_peaks.rds
#
# - Analysis results:
#   * results/tables/overlap_statistics.csv
#   * results/tables/overlap_summary.csv
#   * results/analysis_summary.txt
#   * results/qc_report.txt
#
# - Visualizations:
#   * results/plots/v5_narrow_h2a_overlap_*.pdf
#   * results/plots/v5_broad_h2a_overlap_*.pdf

# ANALYSIS STEPS:
# 1. Load and validate input files
# 2. Process and standardize peak coordinates
# 3. Normalize ChIP-seq signals
# 4. Analyze peak overlaps
# 5. Generate QC metrics
# 6. Create visualizations
# 7. Output summary statistics and reports

# Load required libraries for genomic analysis, visualization and data processing
suppressPackageStartupMessages({
    # Genomic data handling packages
    library(GenomicRanges)    # For handling genomic intervals and operations
    library(rtracklayer)      # For importing/exporting genomic files (BED, bigWig)
    library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome sequence and annotations
    library(GenomeInfoDb)     # Genome information and chromosome standardization
    
    # Peak analysis packages
    library(ChIPseeker)       # ChIP peak annotation and visualization
    library(DiffBind)         # Differential binding analysis for ChIP-seq
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Gene annotations for hg38
    library(org.Hs.eg.db)     # Gene ID mappings for human genes
    
    # Visualization packages
    library(ggplot2)          # For creating publication-quality plots
    library(VennDiagram)      # For visualizing overlaps between peak sets
    library(RColorBrewer)     # For color palettes in visualizations
    library(pheatmap)         # For creating heatmaps
    library(EnrichedHeatmap)  # For ChIP signal profile heatmaps
    library(gridExtra)        # For arranging multiple plots
    
    # Data manipulation packages
    library(dplyr)            # For efficient data manipulation
    library(clusterProfiler)  # For functional enrichment analysis
})

# Define file paths and directory structure
# Base directories for different datasets and their organization
h2a_base_dir <- "../SRF_H2AK119Ub/1_iterative_processing/analysis"
h2a_peaks_dir <- file.path(h2a_base_dir, "peaks")
v5_narrow_peaks_file <- "../SRF_V5/results/peaks/SES-V5ChIP-Seq2_S6_narrow_peaks.narrowPeak"
v5_broad_peaks_file <- "../SRF_V5/results/peaks/SES-V5ChIP-Seq2_S6_broad_peaks.broadPeak"
output_dir <- "results"

# Define peak file paths for H2AK119Ub data
# Organized by peak type (broad), condition (GFP/YAF), and replicate
h2a_peak_files <- list(
    broad = list(
        GFP = list(
            rep1 = file.path(h2a_peaks_dir, "GFP_1_broad_peaks_final.broadPeak"),
            rep2 = file.path(h2a_peaks_dir, "GFP_2_broad_peaks_final.broadPeak"),
            rep3 = file.path(h2a_peaks_dir, "GFP_3_broad_peaks_final.broadPeak")
        ),
        YAF = list(
            rep1 = file.path(h2a_peaks_dir, "YAF_1_broad_peaks_final.broadPeak"),
            rep2 = file.path(h2a_peaks_dir, "YAF_2_broad_peaks_final.broadPeak"),
            rep3 = file.path(h2a_peaks_dir, "YAF_3_broad_peaks_final.broadPeak")
        )
    )
)

# Define bigWig file paths for signal visualization
h2a_bigwig_dir <- file.path(h2a_base_dir, "visualization")
v5_bigwig_dir <- "../SRF_V5/results/bigwig"

# H2AK119Ub bigWig files organized by condition and replicate
h2a_bigwig_files <- list(
    GFP = list(
        rep1 = file.path(h2a_bigwig_dir, "GFP_1.bw"),
        rep2 = file.path(h2a_bigwig_dir, "GFP_2.bw"),
        rep3 = file.path(h2a_bigwig_dir, "GFP_3.bw")
    ),
    YAF = list(
        rep1 = file.path(h2a_bigwig_dir, "YAF_1.bw"),
        rep2 = file.path(h2a_bigwig_dir, "YAF_2.bw"),
        rep3 = file.path(h2a_bigwig_dir, "YAF_3.bw")
    )
)

# V5 ChIP and Input bigWig files
v5_bigwig_files <- list(
    ChIP = file.path(v5_bigwig_dir, "SES-V5ChIP-Seq2_S6.bw"),
    Input = file.path(v5_bigwig_dir, "InputSES-V5ChIP-Seq_S2.bw")
)

# Create output directory structure
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Function to standardize chromosome names and ensure proper genome coordinates
# This is crucial for consistent analysis across different data types
# Function to standardize chromosome names and ensure proper genome coordinates
standardize_chromosomes <- function(gr) {
    # Define standard chromosomes (chr1-22, X, Y)
    std_chroms <- paste0("chr", c(1:22, "X", "Y"))
    
    # First ensure all chromosome names have 'chr' prefix
    current_seqlevels <- seqlevels(gr)
    new_seqlevels <- current_seqlevels
    new_seqlevels <- sub("^chr", "", new_seqlevels)  # Remove any existing 'chr' prefix
    new_seqlevels <- paste0("chr", new_seqlevels)    # Add 'chr' prefix
    
    # Update seqlevels
    seqlevels(gr) <- new_seqlevels
    
    # Get proper sequence info from human genome
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    std_seqinfo <- seqinfo(genome)[std_chroms]
    
    # Keep only standard chromosomes that exist in the input data
    valid_chroms <- intersect(seqlevels(gr), std_chroms)
    if (length(valid_chroms) == 0) {
        stop("No valid chromosomes found after standardization")
    }
    
    # Filter for valid chromosomes
    gr <- keepSeqlevels(gr, valid_chroms, pruning.mode="coarse")
    
    # Update seqinfo
    seqinfo(gr) <- std_seqinfo[seqlevels(gr)]
    
    return(gr)
}

# Function to read and process peak files from multiple replicates
read_peaks <- function(peak_files, type) {
    peaks_list <- list()
    for (condition in names(peak_files[[type]])) {
        condition_peaks <- list()
        for (rep in names(peak_files[[type]][[condition]])) {
            file <- peak_files[[type]][[condition]][[rep]]
            
            # Read the peak file with more flexible numeric parsing
            peaks <- tryCatch({
                # Read the file as a table first
                df <- read.table(file, 
                               colClasses = c("character", "numeric", "numeric", 
                                            "character", "numeric", "character",
                                            "numeric", "numeric", "numeric"),
                               stringsAsFactors = FALSE)
                
                # Convert scientific notation to integers
                df[,2] <- round(as.numeric(df[,2]))
                df[,3] <- round(as.numeric(df[,3]))
                
                # Replace '.' with '*' for strand information
                df[,6] <- ifelse(df[,6] == ".", "*", df[,6])
                
                # Create GRanges object
                GRanges(
                    seqnames = df[,1],
                    ranges = IRanges(start = df[,2], end = df[,3]),
                    strand = df[,6],
                    name = df[,4],
                    score = df[,5],
                    signalValue = df[,7],
                    pValue = df[,8],
                    qValue = df[,9]
                )
            }, error = function(e) {
                warning(sprintf("Error reading file %s: %s", file, e$message))
                return(GRanges())
            })
            
            if (length(peaks) > 0) {
                peaks <- standardize_chromosomes(peaks)
                condition_peaks[[rep]] <- peaks
            }
        }
        peaks_list[[condition]] <- condition_peaks
    }
    return(peaks_list)
}

# Function to ensure genomic coordinates are valid
# Prevents issues with coordinates outside chromosome bounds
ensure_valid_ranges <- function(gr) {
    # Maximum allowed coordinate (2^31 - 1)
    max_coord <- .Machine$integer.max
    
    # Get chromosome lengths from genome info
    chr_lengths <- seqlengths(gr)
    
    # Ensure start positions are between 1 and max allowed
    start(gr) <- pmin(pmax(start(gr), 1), max_coord)
    
    # Ensure end positions are valid and after start positions
    end(gr) <- pmin(pmax(end(gr), start(gr)), max_coord)
    
    # Trim peaks to chromosome lengths if available
    if (!all(is.na(chr_lengths))) {
        gr <- trim(gr)
    }
    
    return(gr)
}

# Function to assess local background signal around peaks
# Helps evaluate peak quality and signal-to-noise ratio
assess_local_background <- function(peaks, bigwig_file, window_size = 10000, min_peaks = 10) {
    # Input validation
    if (length(peaks) < min_peaks) {
        warning("Too few peaks provided for background assessment")
        return(list(
            signal_to_noise = NA,
            qc_metrics = list(
                n_peaks = length(peaks),
                fraction_peaks_with_signal = NA
            )
        ))
    }
    
    # Ensure peaks have valid coordinates
    peaks <- ensure_valid_ranges(peaks)
    
    # Create flanking regions for background assessment
    flank_left <- flank(peaks, width = window_size, start = TRUE)
    flank_right <- flank(peaks, width = window_size, start = FALSE)
    
    # Ensure flanking regions are valid
    flank_left <- ensure_valid_ranges(flank_left)
    flank_right <- ensure_valid_ranges(flank_right)
    
    # Import signal from bigWig file
    signal <- import(bigwig_file, as = "RleList")
    
    # Standardize chromosome naming
    names(signal) <- sub("^chr", "", names(signal))
    names(signal) <- paste0("chr", names(signal))
    
    # Keep only standard chromosomes
    std_chroms <- paste0("chr", c(1:22, "X", "Y"))
    signal <- signal[names(signal) %in% std_chroms]
    peaks <- keepSeqlevels(peaks, std_chroms, pruning.mode="coarse")
    
    # Helper function to safely extract signal values
    safe_extract_signal <- function(regions, signal) {
        tryCatch({
            regions <- trim(regions)
            values <- unlist(Views(signal[seqnames(regions)], 
                                 ranges(regions)))
            values[is.na(values)] <- 0
            return(values)
        }, error = function(e) {
            warning("Error extracting signal: ", e$message)
            return(numeric(0))
        })
    }
    
    # Extract signal values for peaks and flanking regions
    peak_signal <- safe_extract_signal(peaks, signal)
    left_signal <- safe_extract_signal(flank_left, signal)
    right_signal <- safe_extract_signal(flank_right, signal)
    
    # Calculate background signal from flanking regions
    background_signal <- c(left_signal, right_signal)
    
    # Calculate signal-to-noise ratio and other metrics
    signal_to_noise <- if (length(peak_signal) > 0 && length(background_signal) > 0) {
        mean(peak_signal + 1) / mean(background_signal + 1)
    } else {
        NA
    }
    
    fraction_with_signal <- if (length(peak_signal) > 0) {
        mean(peak_signal > 0)
    } else {
        NA
    }
    
    return(list(
        signal_to_noise = signal_to_noise,
        qc_metrics = list(
            n_peaks = length(peaks),
            fraction_peaks_with_signal = fraction_with_signal
        )
    ))
}

# Function to calculate statistical significance of peak overlaps
# Uses permutation testing to assess overlap significance
calculate_overlap_significance <- function(peaks1, peaks2, genome_size = 3.2e9, 
                                        n_permutations = 1000, min_overlap = 5) {
    # Input validation
    if (length(peaks1) == 0 || length(peaks2) == 0) {
        stop("Empty peak sets provided")
    }
    
    # Ensure peaks have valid coordinates
    peaks1 <- ensure_valid_ranges(peaks1)
    peaks2 <- ensure_valid_ranges(peaks2)
    
    # Calculate observed overlap between peak sets
    overlaps <- findOverlaps(peaks1, peaks2)
    observed <- length(unique(queryHits(overlaps)))
    
    if (observed < min_overlap) {
        warning(sprintf("Very few overlaps found (%d < %d)", observed, min_overlap))
    }
    
    # Calculate effective genome size for permutations
    effective_size <- min(genome_size, sum(as.numeric(width(peaks1))))
    
    # Perform permutation test
    permuted_overlaps <- numeric(n_permutations)
    
    for (i in 1:n_permutations) {
        # Generate random regions matching peaks2 properties
        random_starts <- round(runif(length(peaks2), 1, effective_size - max(width(peaks2))))
        random_peaks <- GRanges(
            seqnames = sample(seqlevels(peaks2), length(peaks2), replace = TRUE),
            ranges = IRanges(
                start = random_starts,
                width = width(peaks2)
            )
        )
        
        # Ensure random peaks have valid coordinates
        random_peaks <- ensure_valid_ranges(random_peaks)
        
        # Calculate overlap with random peaks
        random_overlaps <- findOverlaps(peaks1, random_peaks)
        permuted_overlaps[i] <- length(unique(queryHits(random_overlaps)))
    }
    
    # Calculate statistical metrics
    expected <- mean(permuted_overlaps)
    pvalue <- mean(permuted_overlaps >= observed)
    fold_enrichment <- observed / expected
    
    return(list(
        observed = observed,
        expected = expected,
        fold_enrichment = fold_enrichment,
        pvalue = pvalue,
        permuted_overlaps = permuted_overlaps
    ))
}

# Function to merge broad peaks with sophisticated handling
# Useful for combining replicate peaks and handling broad domains
merge_broad_peaks <- function(peaks, min_gap = 1000, min_length = 1000) {
    # Input validation
    if (length(peaks) == 0) {
        warning("Empty peak set provided")
        return(peaks)
    }
    
    # Sort peaks by genomic position
    peaks <- sort(peaks)
    
    # Merge peaks within specified gap distance
    merged <- reduce(peaks, min.gapwidth = min_gap)
    
    # Filter merged peaks by minimum length
    merged <- merged[width(merged) >= min_length]
    
    # Calculate mean signal for merged peaks if score is present
    if ("score" %in% names(mcols(peaks))) {
        overlaps <- findOverlaps(merged, peaks)
        scores <- tapply(mcols(peaks)$score[subjectHits(overlaps)],
                        queryHits(overlaps),
                        mean)
        mcols(merged)$score <- scores[as.character(seq_along(merged))]
    }
    
    # Add metadata about merging process
    mcols(merged)$n_merged <- countOverlaps(merged, peaks)
    mcols(merged)$merged_width <- width(merged)
    
    return(merged)
}

# Function to process peaks with quality control and filtering
process_peaks <- function(peaks_list, is_broad = FALSE, blacklist = NULL) {
    processed_peaks <- list()
    
    for (condition in names(peaks_list)) {
        # Merge replicates into single peak set
        condition_peaks <- unlist(GRangesList(peaks_list[[condition]]))
        
        # Standardize chromosomes and ensure valid coordinates
        condition_peaks <- standardize_chromosomes(condition_peaks)
        condition_peaks <- ensure_valid_ranges(condition_peaks)
        
        if (is_broad) {
            # Use specialized merging for broad peaks
            condition_peaks <- merge_broad_peaks(condition_peaks,
                                              min_gap = 1000,
                                              min_length = 1000)
        } else {
            # Standard peak merging for narrow peaks
            condition_peaks <- reduce(condition_peaks)
        }
        
        # Filter against blacklist regions if provided
        if (!is.null(blacklist)) {
            condition_peaks <- filter_blacklist_regions(condition_peaks, blacklist)
        }
        
        processed_peaks[[condition]] <- condition_peaks
    }
    
    return(processed_peaks)
}

# Function to download and cache ENCODE blacklist regions
# These regions are known problematic regions in the genome
get_blacklist_regions <- function(genome = "hg38", cache_dir = "cache") {
    # Create cache directory
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, paste0(genome, "_blacklist.bed"))
    
    if (!file.exists(cache_file)) {
        # Download blacklist regions if not cached
        blacklist_url <- "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
        download.file(blacklist_url, destfile = paste0(cache_file, ".gz"))
        system(paste("gunzip", paste0(cache_file, ".gz")))
    }
    
    # Read and standardize blacklist regions
    blacklist <- import(cache_file)
    blacklist <- standardize_chromosomes(blacklist)
    
    return(blacklist)
}

# Function to filter peaks against blacklist regions
filter_blacklist_regions <- function(peaks, blacklist) {
    # Find peaks overlapping blacklist regions
    overlaps <- findOverlaps(peaks, blacklist)
    
    # Remove blacklisted peaks
    filtered_peaks <- peaks[-queryHits(overlaps)]
    
    # Add filtering metadata
    mcols(filtered_peaks)$blacklist_filtered <- TRUE
    
    return(filtered_peaks)
}

# Function to validate input files before processing
validate_inputs <- function(chip_bw, input_bw) {
    # Check file existence
    if (!file.exists(chip_bw)) stop(sprintf("ChIP bigwig file not found: %s", chip_bw))
    if (!file.exists(input_bw)) stop(sprintf("Input bigwig file not found: %s", input_bw))
    
    # Validate file format and content
    tryCatch({
        chip_gr <- import(chip_bw, as = "GRanges")
        input_gr <- import(input_bw, as = "GRanges")
    }, error = function(e) {
        stop("Error reading bigwig files: ", e$message)
    })
    
    # Check for empty files
    if (length(chip_gr) == 0) stop("ChIP bigwig file is empty")
    if (length(input_gr) == 0) stop("Input bigwig file is empty")
    
    return(TRUE)
}

# Function to normalize ChIP-seq signal against input
normalize_bigwig_signal <- function(chip_bw, input_bw, scaling_factor = 1, min_reads = 10) {
    # Validate input files
    validate_inputs(chip_bw, input_bw)
    
    # Import ChIP and input signals
    chip_signal <- import(chip_bw, as="RleList")
    input_signal <- import(input_bw, as="RleList")
    
    # Standardize chromosome names
    std_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
    
    # Add 'chr' prefix if missing
    names(chip_signal) <- ifelse(!grepl("^chr", names(chip_signal)), 
                                paste0("chr", names(chip_signal)), 
                                names(chip_signal))
    names(input_signal) <- ifelse(!grepl("^chr", names(input_signal)), 
                                 paste0("chr", names(input_signal)), 
                                 names(input_signal))
    
    # Keep only standard chromosomes
    chip_signal <- chip_signal[names(chip_signal) %in% std_chromosomes]
    input_signal <- input_signal[names(input_signal) %in% std_chromosomes]
    
    # Ensure signals have matching chromosomes
    common_chroms <- intersect(names(chip_signal), names(input_signal))
    chip_signal <- chip_signal[common_chroms]
    input_signal <- input_signal[common_chroms]
    
    # Handle missing values
    chip_signal[is.na(chip_signal)] <- 0
    input_signal[is.na(input_signal)] <- 0
    
    # Apply minimum read count threshold
    input_signal[input_signal < min_reads] <- min_reads
    
    # Calculate normalized signal
    normalized <- (chip_signal + 1)/(input_signal + 1)
    normalized <- normalized * scaling_factor
    
    return(normalized)
}

# Function to categorize peaks based on overlaps
get_peak_categories <- function(v5_peaks, h2a_peaks, prefix) {
    # Find overlaps between V5 and H2AK119Ub peaks
    overlaps <- findOverlaps(v5_peaks, h2a_peaks)
    
    # V5 peaks with H2AK119Ub
    v5_with_h2a <- v5_peaks[unique(queryHits(overlaps))]
    
    # V5 peaks without H2AK119Ub
    v5_only <- v5_peaks[-unique(queryHits(overlaps))]
    
    # H2AK119Ub peaks with V5
    h2a_with_v5 <- h2a_peaks[unique(subjectHits(overlaps))]
    
    # H2AK119Ub peaks without V5
    h2a_only <- h2a_peaks[-unique(subjectHits(overlaps))]
    
    # Create list of categorized peaks
    categorized_peaks <- list(
        v5_with_h2a = v5_with_h2a,
        v5_only = v5_only,
        h2a_with_v5 = h2a_with_v5,
        h2a_only = h2a_only
    )
    
    # Save categorized peaks
    saveRDS(categorized_peaks, 
            file = file.path(output_dir, sprintf("%s_categorized_peaks.rds", prefix)))
    
    return(categorized_peaks)
}

# Main analysis function coordinating all processing steps
main_analysis <- function(narrow_peaks, broad_peaks, v5_narrow_peaks, v5_broad_peaks) {
    # Set up logging
    log_file <- file.path(output_dir, "analysis.log")
    log_connection <- file(log_file, open = "wt")
    on.exit(close(log_connection))
    
    # Helper function for logging
    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        message <- sprintf("[%s] %s", timestamp, msg)
        cat(message, "\n", file = log_connection, append = TRUE)
        cat(message, "\n")
    }
    
    tryCatch({
        log_message("Starting analysis...")
        
        # Create output directories
        dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        
        # Get blacklist regions
        log_message("Loading blacklist regions...")
        blacklist <- get_blacklist_regions()
        
        # Process broad peaks
        log_message("Processing broad peaks...")
        processed_broad <- process_peaks(broad_peaks, is_broad = TRUE, blacklist = blacklist)
        
        # Normalize V5 signal
        log_message("Normalizing V5 signal...")
        normalized_v5_signal <- normalize_bigwig_signal(
            v5_bigwig_files$ChIP,
            v5_bigwig_files$Input,
            min_reads = 10
        )
        
        # Save processed data
        saveRDS(list(
            v5_narrow = v5_narrow_peaks,
            v5_broad = v5_broad_peaks,
            h2a_broad = processed_broad
        ), file = file.path(output_dir, "processed_peaks.rds"))
        
        # Analyze peak overlaps
        narrow_broad_cats <- get_peak_categories(v5_narrow_peaks, processed_broad$YAF, "narrow_broad")
        broad_broad_cats <- get_peak_categories(v5_broad_peaks, processed_broad$YAF, "broad_broad")
        
        # Save normalized signal
        saveRDS(normalized_v5_signal, file = file.path(output_dir, "normalized_v5_signal.rds"))
        
        # Save overlap statistics
        log_message("Saving overlap statistics...")
        overlap_df <- data.frame(
            Comparison = c(
                "V5_narrow_vs_H2A_broad_GFP", "V5_narrow_vs_H2A_broad_YAF",
                "V5_broad_vs_H2A_broad_GFP", "V5_broad_vs_H2A_broad_YAF"
            ),
            Observed = c(
                length(findOverlaps(v5_narrow_peaks, processed_broad$GFP)),
                length(findOverlaps(v5_narrow_peaks, processed_broad$YAF)),
                length(findOverlaps(v5_broad_peaks, processed_broad$GFP)),
                length(findOverlaps(v5_broad_peaks, processed_broad$YAF))
            )
        )
        write.csv(overlap_df, 
                 file = file.path(output_dir, "tables", "overlap_statistics.csv"),
                 row.names = FALSE)
        
        # Generate visualizations
        log_message("Generating plots...")
        generate_plots(NULL, processed_broad, 
                      v5_narrow_peaks, v5_broad_peaks, 
                      normalized_v5_signal, list())
        
        # Generate QC report
        log_message("Generating QC report...")
        generate_qc_report(NULL, processed_broad, 
                         list(), list(), list(), list(), output_dir)
        
        log_message("Analysis completed successfully")
        
    }, error = function(e) {
        log_message(sprintf("ERROR: %s", e$message))
        stop(e)
    })
}

# Function to generate comprehensive QC report
generate_qc_report <- function(narrow_peaks, broad_peaks, narrow_bg_assessment, broad_bg_assessment, narrow_overlap_stats, broad_overlap_stats, output_dir) {
    sink(file.path(output_dir, "qc_report.txt"))
    cat("Quality Control Report\n")
    cat("====================\n\n")
    
    # Helper function to summarize peak analysis results
    summarize_results <- function(peaks, bg_assessment, overlap_stats, peak_type) {
        cat(sprintf("\n%s Peaks Analysis:\n", peak_type))
        cat(sprintf("%s\n", paste(rep("-", nchar(peak_type) + 15), collapse="")))
        
        for (condition in names(peaks)) {
            cat(sprintf("\n%s condition:\n", condition))
            
            # Report background metrics
            bg <- bg_assessment[[condition]]
            cat(sprintf("Background metrics:\n"))
            cat(sprintf("- Signal-to-noise ratio: %.2f\n", bg$signal_to_noise))
            cat(sprintf("- Number of peaks: %d\n", bg$qc_metrics$n_peaks))
            cat(sprintf("- Peaks with signal: %.1f%%\n", 
                       100 * bg$qc_metrics$fraction_peaks_with_signal))
            
            # Report overlap statistics
            ov <- overlap_stats[[condition]]
            cat(sprintf("\nOverlap with V5:\n"))
            cat(sprintf("- Observed overlaps: %d\n", ov$observed))
            cat(sprintf("- Expected overlaps: %.1f\n", ov$expected))
            cat(sprintf("- Fold enrichment: %.2f\n", ov$fold_enrichment))
            cat(sprintf("- P-value: %.2e\n", ov$pvalue))
        }
    }
    
    # Generate summaries for both narrow and broad peaks
    summarize_results(narrow_peaks, narrow_bg_assessment, narrow_overlap_stats, "Narrow")
    summarize_results(broad_peaks, broad_bg_assessment, broad_overlap_stats, "Broad")
    
    sink()
}

# Function to generate analysis plots
generate_plots <- function(narrow_peaks, broad_peaks, v5_narrow_peaks, v5_broad_peaks, normalized_v5_signal, qc_metrics) {
    tryCatch({
        plots_dir <- file.path(output_dir, "plots")
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Helper function for overlap analysis
        analyze_peak_overlap <- function(h2a_peaks, v5_peaks, peak_type, v5_type) {
            overlaps <- findOverlaps(v5_peaks, h2a_peaks)
            v5_categories <- rep("V5 only", length(v5_peaks))
            v5_categories[unique(queryHits(overlaps))] <- "V5 + H2AK119Ub"
            
            peak_data <- data.frame(
                category = v5_categories,
                width = width(v5_peaks),
                signal = numeric(length(v5_peaks))
            )
            
            # Create overlap barplot
            p <- ggplot(peak_data, aes(x = category)) +
                geom_bar() +
                theme_minimal() +
                labs(title = sprintf("%s V5 Peak Overlap with H2AK119Ub (%s)", v5_type, peak_type),
                     x = "Category",
                     y = "Number of Peaks")
            
            # Save plot
            ggsave(file.path(plots_dir, 
                           sprintf("v5_%s_h2a_overlap_%s.pdf", 
                                 tolower(v5_type), 
                                 tolower(peak_type))), 
                   p, width = 8, height = 6)
            
            return(list(
                total_v5 = length(v5_peaks),
                total_h2a = length(h2a_peaks),
                overlapping = length(unique(queryHits(overlaps))),
                peak_data = peak_data
            ))
        }
        
        # Analyze peak overlaps
        narrow_v5_broad_h2a <- analyze_peak_overlap(broad_peaks$YAF, v5_narrow_peaks, "Narrow", "Narrow")
        broad_v5_broad_h2a <- analyze_peak_overlap(broad_peaks$YAF, v5_broad_peaks, "Broad", "Broad")
        
        # Generate summary statistics
        summary_stats <- data.frame(
            Category = c(
                "Total V5 Narrow Peaks",
                "Total V5 Broad Peaks",
                "Total H2AK119Ub Broad Peaks",
                "Narrow V5 + Broad H2AK119Ub",
                "Broad V5 + Broad H2AK119Ub"
            ),
            Count = c(
                length(v5_narrow_peaks),
                length(v5_broad_peaks),
                length(broad_peaks$YAF),
                narrow_v5_broad_h2a$overlapping,
                broad_v5_broad_h2a$overlapping
            )
        )
        
        # Save summary statistics
        write.csv(summary_stats,
                 file = file.path(output_dir, "tables", "overlap_summary.csv"),
                 row.names = FALSE)
        
        # Generate analysis summary
        summary_text <- c(
            "Cross-reference Analysis Summary",
            "==============================",
            "",
            sprintf("Analysis completed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            "",
            "Peak Statistics:",
            sprintf("Total V5 narrow peaks: %d", length(v5_narrow_peaks)),
            sprintf("Total V5 broad peaks: %d", length(v5_broad_peaks)),
            sprintf("H2AK119Ub broad peaks: %d", length(broad_peaks$YAF)),
            "",
            "Overlap Statistics:",
            sprintf("Narrow V5 + Broad H2AK119Ub: %d (%.1f%%)", 
                    narrow_v5_broad_h2a$overlapping,
                    100 * narrow_v5_broad_h2a$overlapping / length(v5_narrow_peaks)),
            sprintf("Broad V5 + Broad H2AK119Ub: %d (%.1f%%)",
                    broad_v5_broad_h2a$overlapping,
                    100 * broad_v5_broad_h2a$overlapping / length(v5_broad_peaks))
        )
        
        writeLines(summary_text, 
                  con = file.path(output_dir, "analysis_summary.txt"))
        
    }, error = function(e) {
        warning(sprintf("Error generating plots: %s", e$message))
        print(e)
    })
}

# Main execution block
main <- function() {
    # Read V5 peaks (both narrow and broad)
    v5_narrow_peaks <- import(v5_narrow_peaks_file)
    v5_narrow_peaks <- standardize_chromosomes(v5_narrow_peaks)
    
    v5_broad_peaks <- import(v5_broad_peaks_file)
    v5_broad_peaks <- standardize_chromosomes(v5_broad_peaks)
    
    # Read only H2AK119Ub broad peaks
    broad_peaks <- read_peaks(h2a_peak_files, "broad")
    
    # Execute main analysis with only broad peaks for H2A
    main_analysis(NULL, broad_peaks, v5_narrow_peaks, v5_broad_peaks)
}

# Execute the main function
main()
