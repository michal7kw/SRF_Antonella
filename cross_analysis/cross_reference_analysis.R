#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ggplot2)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(DiffBind)
    library(VennDiagram)
    library(RColorBrewer)
    library(dplyr)
    library(gridExtra)
    library(clusterProfiler)
    library(pheatmap)
    library(GenomeInfoDb)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(EnrichedHeatmap)
})

# Set paths
h2a_base_dir <- "../SRF_H2AK119Ub/1_iterative_processing/analysis"
h2a_peaks_dir <- file.path(h2a_base_dir, "peaks")
v5_peaks_file <- "../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_narrow_peaks.narrowPeak"
output_dir <- "results"

# Add peak file paths
h2a_peak_files <- list(
    narrow = list(
        GFP = list(
            rep1 = file.path(h2a_peaks_dir, "GFP_1_narrow_peaks.narrowPeak"),
            rep2 = file.path(h2a_peaks_dir, "GFP_2_narrow_peaks.narrowPeak"),
            rep3 = file.path(h2a_peaks_dir, "GFP_3_narrow_peaks.narrowPeak")
        ),
        YAF = list(
            rep1 = file.path(h2a_peaks_dir, "YAF_1_narrow_peaks.narrowPeak"),
            rep2 = file.path(h2a_peaks_dir, "YAF_2_narrow_peaks.narrowPeak"),
            rep3 = file.path(h2a_peaks_dir, "YAF_3_narrow_peaks.narrowPeak")
        )
    ),
    broad = list(
        GFP = list(
            rep1 = file.path(h2a_peaks_dir, "GFP_1_broad_peaks.broadPeak"),
            rep2 = file.path(h2a_peaks_dir, "GFP_2_broad_peaks.broadPeak"),
            rep3 = file.path(h2a_peaks_dir, "GFP_3_broad_peaks.broadPeak")
        ),
        YAF = list(
            rep1 = file.path(h2a_peaks_dir, "YAF_1_broad_peaks.broadPeak"),
            rep2 = file.path(h2a_peaks_dir, "YAF_2_broad_peaks.broadPeak"),
            rep3 = file.path(h2a_peaks_dir, "YAF_3_broad_peaks.broadPeak")
        )
    )
)

# Add bigwig paths
h2a_bigwig_dir <- file.path(h2a_base_dir, "visualization")
v5_bigwig_dir <- "../SRF_V5/bigwig"

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

v5_bigwig_files <- list(
    ChIP = file.path(v5_bigwig_dir, "SES-V5ChIP-Seq2_S6.bw"),
    Input = file.path(v5_bigwig_dir, "InputSES-V5ChIP-Seq_S2.bw")
)

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Enhanced assess_local_background function
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
    
    # Import signal
    signal <- import(bigwig_file, as = "RleList")
    
    # Ensure chromosome naming consistency
    names(signal) <- sub("^chr", "", names(signal))
    names(signal) <- paste0("chr", names(signal))
    
    # Keep only standard chromosomes
    std_chroms <- paste0("chr", c(1:22, "X", "Y"))
    signal <- signal[names(signal) %in% std_chroms]
    peaks <- keepSeqlevels(peaks, std_chroms, pruning.mode="coarse")
    
    # Create flanking regions
    flank_left <- flank(peaks, width = window_size, start = TRUE)
    flank_right <- flank(peaks, width = window_size, start = FALSE)
    
    # Ensure all regions have proper seqinfo
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    std_seqinfo <- seqinfo(genome)[std_chroms]
    
    seqinfo(peaks) <- std_seqinfo[seqlevels(peaks)]
    seqinfo(flank_left) <- std_seqinfo[seqlevels(flank_left)]
    seqinfo(flank_right) <- std_seqinfo[seqlevels(flank_right)]
    
    # Function to safely extract signal
    safe_extract_signal <- function(regions, signal) {
        tryCatch({
            # Ensure regions are within chromosome bounds
            regions <- trim(regions)
            # Extract signal
            values <- unlist(Views(signal[seqnames(regions)], 
                                 ranges(regions)))
            # Replace NAs with 0
            values[is.na(values)] <- 0
            return(values)
        }, error = function(e) {
            warning("Error extracting signal: ", e$message)
            return(numeric(0))
        })
    }
    
    # Extract signal values
    peak_signal <- safe_extract_signal(peaks, signal)
    left_signal <- safe_extract_signal(flank_left, signal)
    right_signal <- safe_extract_signal(flank_right, signal)
    
    # Calculate background as mean of flanking regions
    background_signal <- c(left_signal, right_signal)
    
    # Calculate metrics
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

# Update standardize_chromosomes function
standardize_chromosomes <- function(gr) {
    # Get standard chromosomes
    std_chroms <- paste0("chr", c(1:22, "X", "Y"))
    
    # Ensure chr prefix
    current_seqlevels <- seqlevels(gr)
    new_seqlevels <- sub("^chr", "", current_seqlevels)
    new_seqlevels <- paste0("chr", new_seqlevels)
    seqlevels(gr) <- new_seqlevels
    
    # Get proper sequence info from BSgenome
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    std_seqinfo <- seqinfo(genome)[std_chroms]
    
    # Keep only standard chromosomes and set proper seqinfo
    gr <- keepSeqlevels(gr, std_chroms, pruning.mode="coarse")
    seqinfo(gr) <- std_seqinfo[seqlevels(gr)]
    
    return(gr)
}

# Function to read and process peak files
read_peaks <- function(peak_files, type) {
    peaks_list <- list()
    for (condition in names(peak_files[[type]])) {
        condition_peaks <- list()
        for (rep in names(peak_files[[type]][[condition]])) {
            file <- peak_files[[type]][[condition]][[rep]]
            peaks <- import(file)
            peaks <- standardize_chromosomes(peaks)
            condition_peaks[[rep]] <- peaks
        }
        peaks_list[[condition]] <- condition_peaks
    }
    return(peaks_list)
}

# Function to validate input files
validate_inputs <- function(chip_bw, input_bw) {
    # Check if files exist
    if (!file.exists(chip_bw)) stop(sprintf("ChIP bigwig file not found: %s", chip_bw))
    if (!file.exists(input_bw)) stop(sprintf("Input bigwig file not found: %s", input_bw))
    
    # Check file formats
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

# Enhanced normalize_bigwig_signal function with error handling
normalize_bigwig_signal <- function(chip_bw, input_bw, scaling_factor = 1, min_reads = 10) {
    # Validate inputs
    validate_inputs(chip_bw, input_bw)
    
    # Import signals
    chip_signal <- import(chip_bw, as="RleList")
    input_signal <- import(input_bw, as="RleList")
    
    # Standardize chromosome names
    std_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
    
    # Add 'chr' prefix if missing from names
    names(chip_signal) <- ifelse(!grepl("^chr", names(chip_signal)), 
                                paste0("chr", names(chip_signal)), 
                                names(chip_signal))
    names(input_signal) <- ifelse(!grepl("^chr", names(input_signal)), 
                                 paste0("chr", names(input_signal)), 
                                 names(input_signal))
    
    # Keep only standard chromosomes
    chip_signal <- chip_signal[names(chip_signal) %in% std_chromosomes]
    input_signal <- input_signal[names(input_signal) %in% std_chromosomes]
    
    # Ensure both signals have the same chromosomes in the same order
    common_chroms <- intersect(names(chip_signal), names(input_signal))
    chip_signal <- chip_signal[common_chroms]
    input_signal <- input_signal[common_chroms]
    
    # Handle missing values
    chip_signal[is.na(chip_signal)] <- 0
    input_signal[is.na(input_signal)] <- 0
    
    # Normalize with minimum read count threshold
    input_signal[input_signal < min_reads] <- min_reads
    
    # Calculate normalized signal
    normalized <- (chip_signal + 1)/(input_signal + 1)
    normalized <- normalized * scaling_factor
    
    return(normalized)
}

# Function to download and cache blacklist regions
get_blacklist_regions <- function(genome = "hg38", cache_dir = "cache") {
    # Create cache directory if it doesn't exist
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, paste0(genome, "_blacklist.bed"))
    
    if (!file.exists(cache_file)) {
        # URL for ENCODE blacklist regions
        blacklist_url <- "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
        
        # Download and cache the file
        download.file(blacklist_url, destfile = paste0(cache_file, ".gz"))
        system(paste("gunzip", paste0(cache_file, ".gz")))
    }
    
    # Read blacklist regions
    blacklist <- import(cache_file)
    
    # Standardize chromosomes
    blacklist <- standardize_chromosomes(blacklist)
    
    return(blacklist)
}

# Function to filter peaks against blacklist regions
filter_blacklist_regions <- function(peaks, blacklist) {
    # Find peaks that overlap with blacklist regions
    overlaps <- findOverlaps(peaks, blacklist)
    
    # Keep only peaks that don't overlap with blacklist
    filtered_peaks <- peaks[-queryHits(overlaps)]
    
    # Add metadata about filtering
    mcols(filtered_peaks)$blacklist_filtered <- TRUE
    
    return(filtered_peaks)
}

# Enhanced peak merging for broad peaks
merge_broad_peaks <- function(peaks, min_gap = 1000, min_length = 1000) {
    # Sort peaks by chromosome and start position
    peaks <- sort(peaks)
    
    # Merge peaks that are within min_gap distance
    merged <- reduce(peaks, min.gapwidth = min_gap)
    
    # Filter by minimum length for broad peaks
    merged <- merged[width(merged) >= min_length]
    
    # Calculate mean signal for merged peaks if score is present
    if ("score" %in% names(mcols(peaks))) {
        # For each merged peak, find overlapping original peaks
        overlaps <- findOverlaps(merged, peaks)
        
        # Calculate mean score for each merged peak
        scores <- tapply(mcols(peaks)$score[subjectHits(overlaps)],
                        queryHits(overlaps),
                        mean)
        
        mcols(merged)$score <- scores[as.character(seq_along(merged))]
    }
    
    return(merged)
}

# Enhanced peak processing function
process_peaks <- function(peaks_list, is_broad = FALSE, blacklist = NULL) {
    processed_peaks <- list()
    
    for (condition in names(peaks_list)) {
        # Merge replicates
        condition_peaks <- unlist(GRangesList(peaks_list[[condition]]))
        
        # Standardize chromosomes
        condition_peaks <- standardize_chromosomes(condition_peaks)
        
        if (is_broad) {
            # Use specialized broad peak merging
            condition_peaks <- merge_broad_peaks(condition_peaks,
                                              min_gap = 1000,  # 1kb gap
                                              min_length = 1000) # 1kb min length
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

# Enhanced calculate_overlap_significance function
calculate_overlap_significance <- function(peaks1, peaks2, genome_size = 3.2e9, 
                                        n_permutations = 1000, min_overlap = 5) {
    # Input validation
    if (length(peaks1) == 0 || length(peaks2) == 0) {
        stop("Empty peak sets provided")
    }
    
    # Calculate observed overlap
    observed_overlap <- length(findOverlaps(peaks1, peaks2))
    
    if (observed_overlap < min_overlap) {
        warning(sprintf("Very few overlaps found (%d < %d)", observed_overlap, min_overlap))
    }
    
    # Generate random permutations with progress tracking
    null_distribution <- numeric(n_permutations)
    
    # Calculate effective genome size (excluding blacklisted regions if available)
    effective_genome_size <- genome_size - sum(width(peaks1)) - sum(width(peaks2))
    
    for(i in 1:n_permutations) {
        if (i %% 100 == 0) message(sprintf("Permutation progress: %d%%", round(i/n_permutations * 100)))
        
        # Randomly shuffle peaks while maintaining width distribution
        random_peaks <- peaks2
        max_start <- effective_genome_size - max(width(peaks2))
        if (max_start <= 0) stop("Peaks too large for effective genome size")
        
        random_starts <- sample(1:max_start, length(peaks2))
        ranges(random_peaks) <- IRanges(start = random_starts, width = width(peaks2))
        null_distribution[i] <- length(findOverlaps(peaks1, random_peaks))
    }
    
    # Calculate statistics
    p_value <- sum(null_distribution >= observed_overlap) / n_permutations
    expected_overlap <- mean(null_distribution)
    fold_enrichment <- observed_overlap / expected_overlap
    
    # Calculate confidence intervals
    ci <- quantile(null_distribution, c(0.025, 0.975))
    
    return(list(
        observed = observed_overlap,
        expected = expected_overlap,
        fold_enrichment = fold_enrichment,
        p_value = p_value,
        null_distribution = null_distribution,
        confidence_interval = ci,
        qc_metrics = list(
            n_peaks1 = length(peaks1),
            n_peaks2 = length(peaks2),
            mean_peak_width1 = mean(width(peaks1)),
            mean_peak_width2 = mean(width(peaks2)),
            effective_genome_size = effective_genome_size
        )
    ))
}

# Enhanced main analysis function with proper error handling and logging
main_analysis <- function() {
    log_file <- file.path(output_dir, "analysis.log")
    log_connection <- file(log_file, open = "wt")
    log_message <- function(msg) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        message <- sprintf("[%s] %s", timestamp, msg)
        cat(message, "\n", file = log_connection, append = TRUE)
        cat(message, "\n")
    }
    
    tryCatch({
        log_message("Starting analysis...")
        
        # Download and load blacklist regions
        log_message("Loading blacklist regions...")
        blacklist <- get_blacklist_regions()
        
        # Create output directories
        dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        
        # Process narrow peaks
        log_message("Processing narrow peaks...")
        processed_narrow <- process_peaks(narrow_peaks, 
                                       is_broad = FALSE, 
                                       blacklist = blacklist)
        
        # Process broad peaks with specialized merging
        log_message("Processing broad peaks...")
        processed_broad <- process_peaks(broad_peaks, 
                                      is_broad = TRUE, 
                                      blacklist = blacklist)
        
        # Normalize V5 signal
        log_message("Normalizing V5 signal...")
        normalized_v5_signal <- normalize_bigwig_signal(
            v5_bigwig_files$ChIP, 
            v5_bigwig_files$Input,
            min_reads = 10
        )
        
        # Background assessment
        log_message("Assessing local background...")
        narrow_bg_assessment <- list()
        broad_bg_assessment <- list()
        for (condition in names(processed_narrow)) {
            narrow_bg_assessment[[condition]] <- assess_local_background(
                processed_narrow[[condition]],
                h2a_bigwig_files[[condition]][[1]],
                window_size = 10000,
                min_peaks = 10
            )
            broad_bg_assessment[[condition]] <- assess_local_background(
                processed_broad[[condition]],
                h2a_bigwig_files[[condition]][[1]],
                window_size = 10000,
                min_peaks = 10
            )
        }
        
        # Overlap analysis
        log_message("Calculating overlap significance...")
        narrow_overlap_stats <- list()
        broad_overlap_stats <- list()
        for (condition in names(processed_narrow)) {
            narrow_overlap_stats[[condition]] <- calculate_overlap_significance(
                processed_narrow[[condition]],
                v5_peaks,
                n_permutations = 1000,
                min_overlap = 5
            )
            broad_overlap_stats[[condition]] <- calculate_overlap_significance(
                processed_broad[[condition]],
                v5_peaks,
                n_permutations = 1000,
                min_overlap = 5
            )
        }
        
        # Generate QC report
        log_message("Generating QC report...")
        generate_qc_report(processed_narrow, processed_broad, narrow_bg_assessment, broad_bg_assessment, narrow_overlap_stats, broad_overlap_stats, output_dir)
        
        # Save results
        log_message("Saving analysis results...")
        saveRDS(
            list(
                narrow = list(
                    peaks = processed_narrow,
                    background = narrow_bg_assessment,
                    overlap = narrow_overlap_stats
                ),
                broad = list(
                    peaks = processed_broad,
                    background = broad_bg_assessment,
                    overlap = broad_overlap_stats
                ),
                metadata = list(
                    analysis_date = Sys.time(),
                    genome_build = "hg38",
                    parameters = list(
                        window_size = 10000,
                        min_peaks = 10,
                        min_overlap = 5,
                        n_permutations = 1000
                    )
                )
            ),
            file = file.path(output_dir, "analysis_results.rds")
        )
        
        log_message("Analysis completed successfully")
        
    }, error = function(e) {
        log_message(sprintf("ERROR: %s", e$message))
        stop(e)
    }, finally = {
        close(log_connection)
    })
}

# Function to generate QC report
generate_qc_report <- function(narrow_peaks, broad_peaks, narrow_bg_assessment, broad_bg_assessment, narrow_overlap_stats, broad_overlap_stats, output_dir) {
    sink(file.path(output_dir, "qc_report.txt"))
    cat("Quality Control Report\n")
    cat("====================\n\n")
    
    # Function to summarize results
    summarize_results <- function(peaks, bg_assessment, overlap_stats, peak_type) {
        cat(sprintf("\n%s Peaks Analysis:\n", peak_type))
        cat(sprintf("%s\n", paste(rep("-", nchar(peak_type) + 15), collapse="")))
        
        for (condition in names(peaks)) {
            cat(sprintf("\n%s condition:\n", condition))
            
            # Background metrics
            bg <- bg_assessment[[condition]]
            cat(sprintf("Background metrics:\n"))
            cat(sprintf("- Signal-to-noise ratio: %.2f\n", bg$signal_to_noise))
            cat(sprintf("- Number of peaks: %d\n", bg$qc_metrics$n_peaks))
            cat(sprintf("- Peaks with signal: %.1f%%\n", 
                       100 * bg$qc_metrics$fraction_peaks_with_signal))
            
            # Overlap statistics
            ov <- overlap_stats[[condition]]
            cat(sprintf("\nOverlap with V5:\n"))
            cat(sprintf("- Observed overlaps: %d\n", ov$observed))
            cat(sprintf("- Expected overlaps: %.1f\n", ov$expected))
            cat(sprintf("- Fold enrichment: %.2f\n", ov$fold_enrichment))
            cat(sprintf("- P-value: %.2e\n", ov$p_value))
        }
    }
    
    summarize_results(narrow_peaks, narrow_bg_assessment, narrow_overlap_stats, "Narrow")
    summarize_results(broad_peaks, broad_bg_assessment, broad_overlap_stats, "Broad")
    
    sink()
}

# Read V5 peaks
v5_peaks <- import(v5_peaks_file)
v5_peaks <- standardize_chromosomes(v5_peaks)

# Read H2AK119Ub peaks
narrow_peaks <- read_peaks(h2a_peak_files, "narrow")
broad_peaks <- read_peaks(h2a_peak_files, "broad")

# Execute main analysis
main_analysis()
