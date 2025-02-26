#!/usr/bin/env Rscript

#SBATCH --time=24:00:00        # Increase time limit to 24 hours
#SBATCH --mem=64G             # Request 64GB memory
#SBATCH --cpus-per-task=8     # Request 8 CPUs
#SBATCH --output=cross_ref_%j.log  # Output log with job ID

#####################################################################
# Cross-reference Analysis Script for V5 and H2AK119Ub ChIP-seq Data
#####################################################################

# DESCRIPTION:
# This script performs comprehensive cross-reference analysis between 
# V5 ChIP-seq and H2AK119Ub ChIP-seq data to identify overlapping binding sites.

# INPUT FILES:
# - V5 ChIP-seq peaks:
#   * /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi/peaks/GSM6008236_SESV5ChIPSeq1S5_broadPeak.bed
#   * /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi/peaks/GSM6008237_SESV5ChIPSeq2S6_broadPeak.bed
#
# - V5 ChIP-seq signal:
#   * /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi/bigwig/GSM6008236_SESV5ChIPSeq1S5.bw
#   * /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi/bigwig/GSM6008237_SESV5ChIPSeq2S6.bw
#
# - H2AK119Ub peaks (broad):
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/peaks/GFP_[1-3]_broad_peaks_final.broadPeak
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/peaks/YAF_[1-3]_broad_peaks_final.broadPeak
#
# - H2AK119Ub signal:
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/visualization/GFP_[1-3].bw
#   * ../SRF_H2AK119Ub/1_iterative_processing/analysis/visualization/YAF_[1-3].bw

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
    library(scales)           # For scaling axis labels
    library(grid)             # For grid graphics (required for Venn diagrams)
    
    # Data manipulation packages
    library(dplyr)            # For efficient data manipulation
    library(clusterProfiler)  # For functional enrichment analysis
})

# Define base directories
base_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
h2a_dir <- file.path(base_dir, "SRF_H2AK119Ub")
v5_dir <- file.path(base_dir, "SRF_SES_V5")

# Output directory
output_dir <- file.path(base_dir, "cross_analysis", "results")

# H2AK119Ub file paths
h2a_peaks_dir <- file.path(h2a_dir, "1_iterative_processing", "analysis", "5_peak_calling")
h2a_bigwig_dir <- file.path(h2a_dir, "1_iterative_processing", "analysis", "OLD", "visualization")

# V5 file paths
v5_peaks_dir <- file.path(v5_dir, "results_data_from_ncbi", "peaks")
v5_bigwig_dir <- file.path(v5_dir, "results", "bigwig")

# Define peak file paths for H2AK119Ub data
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

# Define peak file paths for V5 data
v5_broad_peaks_files <- list(
    rep1 = file.path(v5_peaks_dir, "GSM6008236_SESV5ChIPSeq1S5_broadPeak.bed"),
    rep2 = file.path(v5_peaks_dir, "GSM6008237_SESV5ChIPSeq2S6_broadPeak.bed")
)

# Define bigWig file paths for H2AK119Ub data
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

# Define bigWig file paths for V5 data
v5_bigwig_files <- list(
    ChIP = list(
        rep1 = file.path(v5_bigwig_dir, "SRR18590296.bw"),
        rep2 = file.path(v5_bigwig_dir, "SRR18590297.bw")
    ),
    Input = file.path(v5_bigwig_dir, "SRR18590303.bw")
)

# Create output directory structure
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Function to standardize chromosome names and ensure proper genome coordinates
# This is crucial for consistent analysis across different data types
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
    # Helper function to read a single peak file
    read_single_peak_file <- function(peak_file) {
        if (!file.exists(peak_file)) {
            stop(sprintf("Peak file not found: %s", peak_file))
        }
        
        # Generate unique identifier prefix from filename
        file_prefix <- sub("[.][^.]*$", "", basename(peak_file))
        
        # Read raw lines first for counting
        raw_lines <- readLines(peak_file)
        log_message(sprintf("Raw number of lines in file: %d", length(raw_lines)), "DEBUG")
        
        peaks <- tryCatch({
            # Parse the file
            peaks_df <- read.table(peak_file, header = FALSE)
            log_message(sprintf("Number of peaks after initial reading: %d", nrow(peaks_df)), "DEBUG")
            
            # Create GRanges
            gr <- GRanges(
                seqnames = peaks_df$V1,
                ranges = IRanges(
                    start = as.integer(peaks_df$V2),
                    end = as.integer(peaks_df$V3)
                ),
                strand = ifelse(peaks_df$V6 == ".", "*", peaks_df$V6)
            )
            
            # Add metadata columns with non-reserved names
            mcols(gr)$peak_name <- paste0(file_prefix, "_peak_", seq_len(nrow(peaks_df)))
            mcols(gr)$peak_score <- as.numeric(peaks_df$V5)
            
            # Standardize chromosomes
            gr_std <- standardize_chromosomes(gr)
            log_message(sprintf("Number of peaks after standardization: %d", length(gr_std)), "DEBUG")
            
            gr_std
        }, error = function(e) {
            log_message(sprintf("Error reading peak file %s: %s", peak_file, e$message), "ERROR")
            stop(e)
        })
        
        return(peaks)
    }
    
    # For single file input (V5 peaks)
    if (is.character(peak_files)) {
        log_message("Reading single peak file...", "DEBUG")
        return(read_single_peak_file(peak_files))
    }
    
    # For list input (H2AK119Ub peaks)
    if (!is.list(peak_files) || !type %in% names(peak_files)) {
        stop("Invalid peak_files input or type not found in peak_files")
    }
    
    peaks_list <- list()
    for (condition in names(peak_files[[type]])) {
        condition_peaks <- list()
        for (rep in names(peak_files[[type]][[condition]])) {
            file <- peak_files[[type]][[condition]][[rep]]
            log_message(sprintf("Processing file for %s %s: %s", condition, rep, file), "DEBUG")
            peaks <- read_single_peak_file(file)
            
            # Add condition and replicate metadata
            mcols(peaks)$condition <- condition
            mcols(peaks)$replicate <- rep
            
            if (length(peaks) > 0) {
                condition_peaks[[rep]] <- peaks
            }
        }
        
        if (length(condition_peaks) > 0) {
            # Convert to GRangesList and merge
            merged_peaks <- unlist(GRangesList(condition_peaks))
            
            # Generate new unique names
            new_names <- paste0(condition, "_peak_", seq_len(length(merged_peaks)))
            names(merged_peaks) <- new_names
            mcols(merged_peaks)$peak_name <- new_names
            
            peaks_list[[condition]] <- merged_peaks
        }
    }
    
    # Return empty GRangesList if no peaks were found
    if (length(peaks_list) == 0) {
        return(GRangesList())
    }
    
    # Convert the list to a GRangesList
    return(GRangesList(peaks_list))
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
merge_broad_peaks <- function(peaks, min_gap = 1000, min_length = 1000) {
    # Input validation
    if (!is(peaks, "GRanges")) {
        stop("Input must be a GRanges object")
    }
    
    log_message(sprintf("Starting peak merging with %d peaks", length(peaks)), "DEBUG")
    log_message(sprintf("Initial width range: %d - %d bp", 
                       min(width(peaks)), 
                       max(width(peaks))), "DEBUG")
    
    # Sort peaks by chromosome and position
    peaks <- sort(peaks)
    
    # Split peaks by chromosome
    peaks_by_chr <- split(peaks, seqnames(peaks))
    merged_peaks_list <- list()
    
    for (chr in names(peaks_by_chr)) {
        chr_peaks <- peaks_by_chr[[chr]]
        log_message(sprintf("Processing chromosome %s with %d peaks", 
                          chr, length(chr_peaks)), "DEBUG")
        
        # Find overlapping or nearby peaks
        nearby <- findOverlaps(chr_peaks, 
                             maxgap = min_gap,
                             ignore.strand = TRUE)
        
        # Create clusters of overlapping peaks
        clusters <- split(subjectHits(nearby), queryHits(nearby))
        
        # Merge peaks within each cluster
        merged_ranges <- GRanges()
        for (cluster in clusters) {
            cluster_peaks <- chr_peaks[unique(cluster)]
            
            # Merge the ranges in the cluster
            merged <- range(cluster_peaks)
            
            # Calculate mean score if available
            if ("peak_score" %in% names(mcols(cluster_peaks))) {
                mcols(merged)$peak_score <- mean(mcols(cluster_peaks)$peak_score)
            }
            
            merged_ranges <- c(merged_ranges, merged)
        }
        
        # Add any singleton peaks (those without overlaps)
        singleton_idx <- setdiff(seq_along(chr_peaks), 
                               unique(queryHits(nearby)))
        if (length(singleton_idx) > 0) {
            merged_ranges <- c(merged_ranges, chr_peaks[singleton_idx])
        }
        
        merged_peaks_list[[chr]] <- merged_ranges
    }
    
    # Combine all chromosomes
    final_peaks <- unlist(GRangesList(merged_peaks_list))
    
    # Filter by minimum length
    final_peaks <- final_peaks[width(final_peaks) >= min_length]
    
    log_message(sprintf("Finished merging. Final peak count: %d", 
                       length(final_peaks)), "DEBUG")
    log_message(sprintf("Final width range: %d - %d bp", 
                       min(width(final_peaks)), 
                       max(width(final_peaks))), "DEBUG")
    
    return(final_peaks)
}

# Function to process peaks with quality control and filtering
process_peaks <- function(peaks_list, is_broad = FALSE, blacklist = NULL) {
    processed_peaks <- list()
    
    for (condition in names(peaks_list)) {
        # Process each replicate with unique identifiers
        log_message(sprintf("Processing condition: %s", condition), "DEBUG")
        log_message(sprintf("Replicate names: %s", paste(names(peaks_list[[condition]]), collapse=", ")), "DEBUG")
        
        # First read all peaks
        all_peaks <- list()
        total_peaks <- 0
        
        for (rep_name in names(peaks_list[[condition]])) {
            log_message(sprintf("Processing replicate: %s", rep_name), "DEBUG")
            peaks <- read_peaks(peaks_list[[condition]][[rep_name]], "broad")
            
            # Debug peak object
            log_message(sprintf("Number of peaks for %s: %d", rep_name, length(peaks)), "DEBUG")
            
            # Add unique identifiers
            mcols(peaks)$peak_id <- paste0(condition, "_", rep_name, "_", seq_len(length(peaks)))
            mcols(peaks)$replicate <- rep_name
            names(peaks) <- mcols(peaks)$peak_id
            
            all_peaks[[rep_name]] <- peaks
            total_peaks <- total_peaks + length(peaks)
        }
        
        # Combine all peaks with unique names
        combined_peaks <- unlist(GRangesList(all_peaks))
        names(combined_peaks) <- paste0(condition, "_peak_", seq_len(length(combined_peaks)))
        mcols(combined_peaks)$peak_id <- names(combined_peaks)
        condition_peaks <- combined_peaks
        
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

# Function to load ENCODE blacklist regions from local file
# These regions are known problematic regions in the genome
get_blacklist_regions <- function(genome = "hg38") {
    # Path to local blacklist file
    blacklist_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/hg38-blacklist.v2.bed"
    
    # Verify file exists
    if (!file.exists(blacklist_file)) {
        stop(sprintf("Blacklist file not found: %s", blacklist_file))
    }
    
    # Read and standardize blacklist regions
    blacklist <- import(blacklist_file)
    blacklist <- standardize_chromosomes(blacklist)
    
    return(blacklist)
}

# Function to filter blacklist regions with more lenient criteria
filter_blacklist_regions <- function(peaks, blacklist) {
    # Find overlaps between peaks and blacklist
    overlaps <- findOverlaps(peaks, blacklist)
    
    # Calculate overlap percentages
    peak_widths <- width(peaks)[queryHits(overlaps)]
    blacklist_widths <- width(blacklist)[subjectHits(overlaps)]
    overlap_widths <- width(pintersect(peaks[queryHits(overlaps)], blacklist[subjectHits(overlaps)]))
    
    # Calculate percentage of peak overlapping with blacklist
    overlap_percent <- overlap_widths / peak_widths
    
    # Only remove peaks that overlap significantly with blacklist (>50% overlap)
    peaks_to_remove <- unique(queryHits(overlaps)[overlap_percent > 0.5])
    
    if (length(peaks_to_remove) > 0) {
        filtered_peaks <- peaks[-peaks_to_remove]
    } else {
        filtered_peaks <- peaks
    }
    
    log_message(sprintf("Filtered %d peaks (%.1f%%) that significantly overlap with blacklist regions", 
                       length(peaks_to_remove), 
                       100 * length(peaks_to_remove) / length(peaks)), 
                "DEBUG")
    
    return(filtered_peaks)
}

# Function to normalize ChIP-seq signal against input
normalize_bigwig_signal <- function(chip_bw_list, input_bw, scaling_factor = 1, min_reads = 10) {
    # Validate input files
    if (!file.exists(input_bw)) stop(sprintf("Input bigwig file not found: %s", input_bw))
    
    # Import input signal
    input_signal <- import(input_bw, as="RleList")
    
    # Process each replicate
    normalized_signals <- lapply(chip_bw_list, function(chip_bw) {
        if (!file.exists(chip_bw)) stop(sprintf("ChIP bigwig file not found: %s", chip_bw))
        
        # Import ChIP signal
        chip_signal <- import(chip_bw, as="RleList")
        
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
    })
    
    # Average the replicates
    averaged_signal <- Reduce("+", normalized_signals) / length(normalized_signals)
    
    return(averaged_signal)
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

# Update the validate_inputs function - remove the malformed line
validate_inputs <- function(chip_bw_list, input_bw) {
    # Check input file
    if (!file.exists(input_bw)) stop(sprintf("Input bigwig file not found: %s", input_bw))
    
    # Check each ChIP replicate file
    for (chip_bw in chip_bw_list) {
        if (!file.exists(chip_bw)) stop(sprintf("ChIP bigwig file not found: %s", chip_bw))
    }
    
    return(TRUE)
}

# Enhanced log_message function with debug levels
log_message <- function(msg, level = "INFO", log_connection = NULL) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message <- sprintf("[%s] [%s] %s", timestamp, level, msg)
    
    # Write to log file if connection is provided
    if (!is.null(log_connection)) {
        cat(message, "\n", file = log_connection, append = TRUE)
    }
    
    # Also print to console
    cat(message, "\n")
}

# Function to validate peak characteristics
validate_peaks <- function(peaks, peak_type = "unknown") {
    # Basic validation
    n_peaks <- length(peaks)
    peak_widths <- width(peaks)
    
    # Calculate key metrics
    avg_width <- mean(peak_widths)
    median_width <- median(peak_widths)
    min_width <- min(peak_widths)
    max_width <- max(peak_widths)
    
    # Chromosome distribution
    chrom_dist <- table(seqnames(peaks))
    
    # Biological validation thresholds
    MIN_EXPECTED_PEAKS <- 1000  # Adjust based on your expectations
    MAX_EXPECTED_PEAKS <- 100000
    MIN_PEAK_WIDTH <- 100
    MAX_PEAK_WIDTH <- 1000000
    
    # Log validation results
    log_message(sprintf("\nValidating %s peaks:", peak_type), "DEBUG")
    log_message(sprintf("Number of peaks: %d", n_peaks), "DEBUG")
    log_message(sprintf("Peak width statistics:"), "DEBUG")
    log_message(sprintf("  - Average: %.2f bp", avg_width), "DEBUG")
    log_message(sprintf("  - Median: %.2f bp", median_width), "DEBUG")
    log_message(sprintf("  - Range: %d - %d bp", min_width, max_width), "DEBUG")
    
    # Log chromosome distribution
    log_message("Chromosome distribution:", "DEBUG")
    for (chrom in names(chrom_dist)) {
        log_message(sprintf("  - %s: %d peaks", chrom, chrom_dist[chrom]), "DEBUG")
    }
    
    # Warnings for potential issues
    if (n_peaks < MIN_EXPECTED_PEAKS) {
        log_message(sprintf("WARNING: Low number of peaks (%d < %d)", n_peaks, MIN_EXPECTED_PEAKS), "WARN")
    }
    if (n_peaks > MAX_EXPECTED_PEAKS) {
        log_message(sprintf("WARNING: Unusually high number of peaks (%d > %d)", n_peaks, MAX_EXPECTED_PEAKS), "WARN")
    }
    if (min_width < MIN_PEAK_WIDTH) {
        log_message(sprintf("WARNING: Some peaks are unusually short (%d bp < %d bp)", min_width, MIN_PEAK_WIDTH), "WARN")
    }
    if (max_width > MAX_PEAK_WIDTH) {
        log_message(sprintf("WARNING: Some peaks are unusually long (%d bp > %d bp)", max_width, MAX_PEAK_WIDTH), "WARN")
    }
    
    return(list(
        n_peaks = n_peaks,
        avg_width = avg_width,
        median_width = median_width,
        width_range = c(min_width, max_width),
        chrom_dist = chrom_dist
    ))
}

# Function to check signal distribution
check_signal_distribution <- function(signal, name = "unknown") {
    # Convert signal to numeric vector
    signal_values <- as.numeric(unlist(signal))
    
    # Calculate statistics
    stats <- list(
        mean = mean(signal_values, na.rm = TRUE),
        median = median(signal_values, na.rm = TRUE),
        sd = sd(signal_values, na.rm = TRUE),
        quantiles = quantile(signal_values, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    )
    
    # Log signal statistics
    log_message(sprintf("\nSignal distribution for %s:", name), "DEBUG")
    log_message(sprintf("Mean signal: %.2f", stats$mean), "DEBUG")
    log_message(sprintf("Median signal: %.2f", stats$median), "DEBUG")
    log_message(sprintf("Standard deviation: %.2f", stats$sd), "DEBUG")
    log_message("Quantiles:", "DEBUG")
    log_message(sprintf("  - 0%%: %.2f", stats$quantiles[1]), "DEBUG")
    log_message(sprintf("  - 25%%: %.2f", stats$quantiles[2]), "DEBUG")
    log_message(sprintf("  - 50%%: %.2f", stats$quantiles[3]), "DEBUG")
    log_message(sprintf("  - 75%%: %.2f", stats$quantiles[4]), "DEBUG")
    log_message(sprintf("  - 100%%: %.2f", stats$quantiles[5]), "DEBUG")
    
    # Check for potential issues
    if (any(is.na(signal_values))) {
        log_message(sprintf("WARNING: %d NA values found in signal", sum(is.na(signal_values))), "WARN")
    }
    if (any(signal_values < 0, na.rm = TRUE)) {
        log_message("WARNING: Negative signal values detected", "WARN")
    }
    
    return(stats)
}

# Function to generate peak statistics summary
generate_peak_stats <- function(peaks, name) {
    # Calculate basic statistics
    widths <- width(peaks)
    stats <- list(
        count = length(peaks),
        mean_width = mean(widths),
        median_width = median(widths),
        min_width = min(widths),
        max_width = max(widths),
        sd_width = sd(widths)
    )
    
    # Calculate chromosome distribution
    chrom_dist <- table(seqnames(peaks))
    
    # Create summary text
    summary <- sprintf("\n=== %s Peak Statistics ===\n", name)
    summary <- paste0(summary, sprintf("Total Peaks: %d\n", stats$count))
    summary <- paste0(summary, "\nPeak Width Statistics:\n")
    summary <- paste0(summary, sprintf("- Mean: %.1f bp\n", stats$mean_width))
    summary <- paste0(summary, sprintf("- Median: %.1f bp\n", stats$median_width))
    summary <- paste0(summary, sprintf("- Range: %d - %d bp\n", stats$min_width, stats$max_width))
    summary <- paste0(summary, sprintf("- Standard Deviation: %.1f bp\n", stats$sd_width))
    
    summary <- paste0(summary, "\nChromosome Distribution:\n")
    for (chr in names(chrom_dist)) {
        summary <- paste0(summary, sprintf("- %s: %d peaks (%.1f%%)\n", 
                                         chr, chrom_dist[chr], 
                                         100 * chrom_dist[chr]/stats$count))
    }
    
    return(list(text = summary, stats = stats, chrom_dist = chrom_dist))
}

# Function to generate filtering summary
generate_filtering_summary <- function(original, filtered, step_name) {
    removed <- length(original) - length(filtered)
    percent_removed <- 100 * removed / length(original)
    
    summary <- sprintf("\n=== %s Filtering Summary ===\n", step_name)
    summary <- paste0(summary, sprintf("Original peaks: %d\n", length(original)))
    summary <- paste0(summary, sprintf("Filtered peaks: %d\n", length(filtered)))
    summary <- paste0(summary, sprintf("Removed peaks: %d (%.1f%%)\n", removed, percent_removed))
    
    if (length(filtered) > 0) {
        width_changes <- list(
            original = width(original),
            filtered = width(filtered)
        )
        
        summary <- paste0(summary, "\nWidth Statistics Changes:\n")
        summary <- paste0(summary, sprintf("Original - Mean: %.1f, Median: %.1f, Range: %d-%d\n",
                                         mean(width_changes$original),
                                         median(width_changes$original),
                                         min(width_changes$original),
                                         max(width_changes$original)))
        summary <- paste0(summary, sprintf("Filtered - Mean: %.1f, Median: %.1f, Range: %d-%d\n",
                                         mean(width_changes$filtered),
                                         median(width_changes$filtered),
                                         min(width_changes$filtered),
                                         max(width_changes$filtered)))
    }
    
    return(summary)
}

# Function to generate overlap analysis report
generate_overlap_report <- function(v5_peaks, h2a_peaks, condition) {
    overlaps <- findOverlaps(v5_peaks, h2a_peaks)
    v5_with_h2a <- unique(queryHits(overlaps))
    h2a_with_v5 <- unique(subjectHits(overlaps))
    
    # Calculate overlap widths
    overlap_widths <- width(pintersect(v5_peaks[queryHits(overlaps)], 
                                     h2a_peaks[subjectHits(overlaps)]))
    
    report <- sprintf("\n=== Overlap Analysis: V5 vs H2AK119Ub %s ===\n", condition)
    report <- paste0(report, sprintf("\nOverall Statistics:\n"))
    report <- paste0(report, sprintf("- Total V5 peaks: %d\n", length(v5_peaks)))
    report <- paste0(report, sprintf("- Total H2AK119Ub peaks: %d\n", length(h2a_peaks)))
    report <- paste0(report, sprintf("- V5 peaks with overlap: %d (%.1f%%)\n",
                                   length(v5_with_h2a),
                                   100 * length(v5_with_h2a)/length(v5_peaks)))
    report <- paste0(report, sprintf("- H2AK119Ub peaks with overlap: %d (%.2f%%)\n",
                                   length(h2a_with_v5),
                                   100 * length(h2a_with_v5)/length(h2a_peaks)))
    
    if (length(overlap_widths) > 0) {
        report <- paste0(report, "\nOverlap Width Statistics:\n")
        report <- paste0(report, sprintf("- Mean: %.1f bp\n", mean(overlap_widths)))
        report <- paste0(report, sprintf("- Median: %.1f bp\n", median(overlap_widths)))
        report <- paste0(report, sprintf("- Range: %d - %d bp\n",
                                       min(overlap_widths),
                                       max(overlap_widths)))
    }
    
    return(report)
}

# Add required libraries for plotting
suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(RColorBrewer)
    library(scales)
})

# Function to plot peak width distributions
plot_peak_width_distribution <- function(peaks_list, names, title = "Peak Width Distribution") {
    # Convert GRanges objects to data frames with widths
    plot_data <- data.frame()
    for (i in seq_along(peaks_list)) {
        # Extract widths from GRanges object
        widths_df <- data.frame(
            width = width(peaks_list[[i]]),
            dataset = names[i]
        )
        plot_data <- rbind(plot_data, widths_df)
    }
    
    # Create density plot
    p <- ggplot(plot_data, aes(x = width, fill = dataset)) +
        geom_density(alpha = 0.5) +
        scale_x_log10(labels = comma) +
        theme_bw() +
        labs(x = "Peak Width (bp)", y = "Density", title = title) +
        theme(legend.position = "bottom")
    
    return(p)
}

# Function to plot chromosome distribution
plot_chromosome_distribution <- function(peaks_list, names, title = "Peak Distribution by Chromosome") {
    # Combine chromosome counts from all peak sets
    plot_data <- do.call(rbind, lapply(seq_along(peaks_list), function(i) {
        chrom_counts <- table(seqnames(peaks_list[[i]]))
        data.frame(
            chromosome = names(chrom_counts),
            count = as.numeric(chrom_counts),
            dataset = names[i]
        )
    }))
    
    # Order chromosomes naturally
    plot_data$chromosome <- factor(plot_data$chromosome, 
                                 levels = paste0("chr", c(1:22, "X", "Y")))
    
    # Create bar plot
    p <- ggplot(plot_data, aes(x = chromosome, y = count, fill = dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        labs(x = "Chromosome", y = "Number of Peaks", title = title)
    
    return(p)
}

# Function to analyze peak overlaps and generate comprehensive statistics
analyze_overlaps <- function(v5_peaks, h2a_peaks, condition) {
    # Find overlaps
    overlaps <- findOverlaps(v5_peaks, h2a_peaks)
    v5_with_h2a <- unique(queryHits(overlaps))
    h2a_with_v5 <- unique(subjectHits(overlaps))
    
    # Calculate overlap statistics
    overlap_stats <- list(
        total_v5 = length(v5_peaks),
        total_h2a = length(h2a_peaks),
        overlapping_peaks = length(v5_with_h2a),
        percent_v5_overlap = 100 * length(v5_with_h2a) / length(v5_peaks),
        percent_h2a_overlap = 100 * length(h2a_with_v5) / length(h2a_peaks)
    )
    
    # Calculate overlap widths
    if (length(overlaps) > 0) {
        overlap_widths <- width(pintersect(v5_peaks[queryHits(overlaps)], 
                                         h2a_peaks[subjectHits(overlaps)]))
        overlap_stats$width_stats <- list(
            min = min(overlap_widths),
            max = max(overlap_widths),
            mean = mean(overlap_widths),
            median = median(overlap_widths),
            sd = sd(overlap_widths)
        )
    }
    
    # Create summary plot focusing on overlapping peaks
    summary_data <- data.frame(
        category = c("Total Peaks", "Overlapping Peaks"),
        V5 = c(length(v5_peaks), length(v5_with_h2a)),
        H2AK119Ub = c(length(h2a_peaks), length(h2a_with_v5))
    )
    
    # Melt data for plotting
    plot_data <- reshape2::melt(summary_data, id.vars = "category")
    
    p <- ggplot(plot_data, aes(x = category, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_log10(labels = comma) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "", y = "Number of Peaks",
             title = sprintf("Peak Analysis Summary - %s", condition),
             fill = "Peak Type")
    
    return(list(stats = overlap_stats, plot = p))
}

# Function to generate and analyze peak profiles with improved methodology
generate_peak_profiles <- function(peaks, signal_files, condition, window_size = 2000, bin_size = 50, 
                                      max_peaks = 10000, chunk_size = 1000) {
    if (length(peaks) == 0) {
        log_message("No peaks provided for profile generation", "WARNING")
        return(NULL)
    }
    
    # Ensure peaks is a GRanges object
    if (!is(peaks, "GRanges")) {
        stop("Input 'peaks' must be a GRanges object")
    }
    
    # Subsample peaks if there are too many
    if (length(peaks) > max_peaks) {
        set.seed(42)  # For reproducibility
        idx <- sample(length(peaks), max_peaks)
        peaks <- peaks[idx]
        log_message(sprintf("Subsampled peaks to %d for profile generation", max_peaks), "INFO")
    }
    
    log_message(sprintf("Generating profiles for %d peaks in %s condition", 
                       length(peaks), condition), "INFO")
    
    # Center peaks and create windows
    peak_windows <- resize(peaks, width = window_size * 2, fix = "center")
    
    # Process each signal file
    profile_data <- list()
    for (signal_file in signal_files) {
        if (!file.exists(signal_file)) {
            log_message(sprintf("Signal file not found: %s", signal_file), "WARNING")
            next
        }
        
        tryCatch({
            # Import and process signal
            signal <- import.bw(signal_file)
            signal <- standardize_chromosomes(signal)
            
            # Calculate coverage matrix using GenomicRanges methods
            coverage_mat <- normalizeToMatrix(signal, peak_windows,
                                            value_column = "score",
                                            mean_mode = "w0",
                                            w = bin_size,
                                            extend = 0,
                                            smooth = TRUE)
            
            if (all(is.na(coverage_mat)) || all(coverage_mat == 0)) {
                log_message(sprintf("No signal detected in %s", basename(signal_file)), "WARNING")
                next
            }
            
            # Calculate summary statistics
            profile_mean <- colMeans(coverage_mat, na.rm = TRUE)
            profile_sd <- apply(coverage_mat, 2, sd, na.rm = TRUE)
            
            profile_data[[basename(signal_file)]] <- list(
                matrix = coverage_mat,
                mean = profile_mean,
                sd = profile_sd
            )
            
            log_message(sprintf("Processed %s successfully", basename(signal_file)), "INFO")
        }, error = function(e) {
            log_message(sprintf("Error processing %s: %s", basename(signal_file), e$message), "ERROR")
        })
    }
    
    if (length(profile_data) == 0) {
        log_message("No valid profiles generated", "WARNING")
        return(NULL)
    }
    
    return(profile_data)
}

# Function to create profile plots from profile data
create_profile_plot <- function(profile_data, condition) {
    # Prepare data for plotting
    plot_data <- lapply(names(profile_data), function(sample) {
        data <- profile_data[[sample]]
        positions <- seq(-1000, 1000, length.out = length(data$mean))
        
        data.frame(
            position = positions,
            signal = data$mean,
            sd_low = data$mean - data$sd,
            sd_high = data$mean + data$sd,
            sample = sample
        )
    })
    
    plot_data <- do.call(rbind, plot_data)
    
    # Create plot
    p <- ggplot(plot_data, aes(x = position, y = signal, color = sample)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = sd_low, ymax = sd_high, fill = sample),
                    alpha = 0.2, color = NA) +
        theme_bw() +
        labs(x = "Distance from Peak Center (bp)",
             y = "Average Signal",
             title = sprintf("%s Peak Profiles", condition)) +
        theme(legend.position = "bottom",
              plot.title = element_text(hjust = 0.5))
    
    return(p)
}

# Main analysis function with enhanced documentation
# Main analysis function incorporating best practices for ChIP-seq analysis
main_analysis <- function(narrow_peaks, broad_peaks, v5_narrow_peaks, v5_broad_peaks) {
    tryCatch({
        # Set up logging and report files
        log_file <- file.path(output_dir, "analysis.log")
        report_file <- file.path(output_dir, "analysis_report.txt")
        log_connection <- file(log_file, open = "wt")
        report_connection <- file(report_file, open = "wt")
        
        # Write analysis header
        cat("ChIP-seq Cross-reference Analysis Report\n",
            "=====================================\n\n",
            sprintf("Analysis Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            file = report_connection)
        
        # Process V5 peaks
        log_message("Processing V5 peaks...", "INFO")
        v5_peaks <- v5_broad_peaks  # Using broad peaks as per best practices
        
        # Write V5 peak statistics
        cat("\nV5 Peak Statistics:\n",
            "------------------\n",
            sprintf("Total peaks: %d\n", length(v5_peaks)),
            sprintf("Width range: %d - %d bp\n", min(width(v5_peaks)), max(width(v5_peaks))),
            sprintf("Median width: %d bp\n", median(width(v5_peaks))),
            file = report_connection)
        
        # Analyze peaks by condition
        for (condition in names(broad_peaks)) {
            log_message(sprintf("Analyzing %s condition", condition), "INFO")
            
            cat(sprintf("\n%s Analysis:\n", condition),
                sprintf("%s\n", paste(rep("-", nchar(condition) + 10), collapse = "")),
                file = report_connection)
            
            # Get H2AK119Ub peaks for this condition
            h2a_peaks <- broad_peaks[[condition]]
            
            # Annotate peaks
            log_message("Performing peak annotation", "INFO")
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            
            # Debug peak object before conversion
            log_message(sprintf("Number of peaks before conversion: %d", length(h2a_peaks)), "DEBUG")
            log_message(sprintf("Names of h2a_peaks: %s", paste(head(names(h2a_peaks)), collapse=", ")), "DEBUG")
            
            # Convert peaks to data frame preserving ALL metadata
            peaks_df <- as.data.frame(h2a_peaks)
            original_metadata <- mcols(h2a_peaks)
            log_message(sprintf("Peaks data frame dimensions: %d x %d", nrow(peaks_df), ncol(peaks_df)), "DEBUG")
            
            # Check for existing metadata
            existing_mcols <- colnames(original_metadata)
            log_message(sprintf("Existing metadata columns: %s", paste(existing_mcols, collapse=", ")), "DEBUG")
            
            # Preserve or generate peak IDs
            if ("peak_id" %in% existing_mcols && 
                !all(is.na(original_metadata$peak_id)) && 
                !all(original_metadata$peak_id == "")) {
                # Use existing peak IDs
                peak_ids <- original_metadata$peak_id
                log_message("Using existing peak IDs from metadata", "DEBUG")
            } else if (!is.null(names(h2a_peaks)) && 
                       !all(is.na(names(h2a_peaks))) && 
                       !all(names(h2a_peaks) == "")) {
                # Use existing names
                peak_ids <- names(h2a_peaks)
                log_message("Using existing peak names", "DEBUG")
            } else {
                # Generate new unique IDs
                peak_ids <- paste0(condition, "_peak_", seq_len(nrow(peaks_df)))
                log_message("Generating new peak IDs", "DEBUG")
            }
            
            # Ensure peak IDs are unique
            if (length(unique(peak_ids)) != length(peak_ids)) {
                log_message("WARNING: Duplicate peak IDs found, generating new unique IDs", "WARNING")
                peak_ids <- make.unique(peak_ids)
            }
            
            # Create GRanges preserving ALL metadata
            peaks_gr <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns=TRUE)
            
            # Copy over ALL original metadata
            if (length(existing_mcols) > 0) {
                for (col in existing_mcols) {
                    mcols(peaks_gr)[[col]] <- original_metadata[[col]]
                }
            }
            
            # Update/add required metadata
            mcols(peaks_gr)$peak_id <- peak_ids
            mcols(peaks_gr)$condition <- condition
            names(peaks_gr) <- peak_ids
            
            # Validate final object
            stopifnot(
                length(peaks_gr) == length(h2a_peaks),
                length(unique(names(peaks_gr))) == length(peaks_gr),
                !any(is.na(names(peaks_gr))),
                !any(names(peaks_gr) == "")
            )
            
            # Debug final object
            log_message(sprintf("Final GRanges metadata columns: %s", 
                paste(colnames(mcols(peaks_gr)), collapse=", ")), "DEBUG")
            
            log_message(sprintf("Final GRanges object length: %d", length(peaks_gr)), "DEBUG")
            log_message(sprintf("Final GRanges names: %s", paste(head(names(peaks_gr)), collapse=", ")), "DEBUG")
            
            # Perform peak annotation
            peak_annot <- annotatePeak(peaks_gr, TxDb=txdb,
                                      annoDb="org.Hs.eg.db")
            
            # Write annotation summary
            annot_summary <- as.data.frame(peak_annot)
            cat(sprintf("\nPeak Annotation Summary:\n"),
                sprintf("Total peaks: %d\n", nrow(annot_summary)),
                sprintf("Promoter peaks: %d (%.1f%%)\n", 
                        sum(grepl("Promoter", annot_summary$annotation)),
                        100 * sum(grepl("Promoter", annot_summary$annotation)) / nrow(annot_summary)),
                file = report_connection)
            
            # Analyze overlaps
            log_message("Analyzing peak overlaps", "INFO")
            overlap_analysis <- analyze_overlaps(v5_peaks, h2a_peaks, condition)
            
            # Write overlap statistics
            cat("\nOverlap Analysis:\n",
                sprintf("V5 peaks overlapping with H2AK119Ub: %d (%.1f%%)\n",
                        overlap_analysis$stats$overlapping_peaks,
                        overlap_analysis$stats$percent_v5_overlap),
                sprintf("H2AK119Ub peaks overlapping with V5: %d (%.1f%%)\n",
                        overlap_analysis$stats$overlapping_peaks,
                        overlap_analysis$stats$percent_h2a_overlap),
                file = report_connection)
            
            if (!is.null(overlap_analysis$stats$width_stats)) {
                cat("\nOverlap Width Statistics:\n",
                    sprintf("Median width: %d bp\n", overlap_analysis$stats$width_stats$median),
                    sprintf("Mean width: %.1f bp\n", overlap_analysis$stats$width_stats$mean),
                    sprintf("Width range: %d - %d bp\n",
                            overlap_analysis$stats$width_stats$min,
                            overlap_analysis$stats$width_stats$max),
                    file = report_connection)
            }
            
            # Generate and save overlap plot
            pdf(file.path(output_dir, "plots", sprintf("%s_overlap_analysis.pdf", condition)))
            print(overlap_analysis$plot)
            dev.off()
            
            # Generate peak profiles
            log_message("Generating peak profiles", "INFO")
            profile_data <- generate_peak_profiles(h2a_peaks, 
                                                 h2a_bigwig_files[[condition]], 
                                                 condition)
            
            if (!is.null(profile_data)) {
                # Create profile plot
                profile_plot <- create_profile_plot(profile_data, condition)
                pdf(file.path(output_dir, "plots", 
                             sprintf("%s_peak_profiles.pdf", condition)))
                print(profile_plot)
                dev.off()
            }
        }
        report_connection <- file(report_file, open = "wt")
        on.exit({
            close(log_connection)
            close(report_connection)
        })
        
        # Write report header
        cat("===========================================\n", file = report_connection)
        cat("ChIP-seq Cross-Reference Analysis Report\n", file = report_connection)
        cat("===========================================\n\n", file = report_connection)
        cat(sprintf("Analysis Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
            file = report_connection)
        
        # Initial V5 peak analysis
        log_message("Analyzing initial V5 peaks...", "INFO", log_connection)
        v5_initial_stats <- generate_peak_stats(v5_broad_peaks, "Initial V5")
        cat(v5_initial_stats$text, file = report_connection)
        
        # Get blacklist regions
        log_message("Loading blacklist regions...", "INFO", log_connection)
        blacklist <- get_blacklist_regions()
        
        # Filter V5 peaks against blacklist
        # log_message("Filtering V5 peaks against blacklist...", "INFO", log_connection)
        # v5_filtered <- filter_blacklist_regions(v5_broad_peaks, blacklist)
        v5_filtered <- v5_broad_peaks
        filtering_summary <- generate_filtering_summary(v5_broad_peaks, v5_filtered, "Blacklist")
        cat(filtering_summary, file = report_connection)
        
        # Analyze filtered V5 peaks
        v5_filtered_stats <- generate_peak_stats(v5_filtered, "Filtered V5")
        cat(v5_filtered_stats$text, file = report_connection)
        
        # Create plots directory
        plots_dir <- file.path(output_dir, "plots")
        dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Plot V5 peak width distribution before and after filtering
        width_dist <- plot_peak_width_distribution(
            list(v5_broad_peaks, v5_filtered),
            c("Original", "Filtered"),
            "V5 Peak Width Distribution"
        )
        ggsave(file.path(plots_dir, "v5_peak_width_distribution.pdf"), width_dist,
               width = 8, height = 6)
        
        # Plot V5 chromosome distribution
        chrom_dist <- plot_chromosome_distribution(
            list(v5_broad_peaks, v5_filtered),
            c("Original", "Filtered"),
            "V5 Peak Distribution by Chromosome"
        )
        ggsave(file.path(plots_dir, "v5_chromosome_distribution.pdf"), chrom_dist,
               width = 10, height = 6)
        
        # Process each H2AK119Ub condition with plots
        for (condition in names(broad_peaks)) {
            cat(sprintf("\n\n=== Analysis for %s condition ===\n", condition), 
                file = report_connection)
            
            # Perform peak annotation and overlap analysis
            log_message(sprintf("Analyzing %s peaks", condition), "INFO")
            
            # Annotate peaks using ChIPseeker
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            peaks_gr <- as(broad_peaks[[condition]], "GRanges")
            peak_annot <- annotatePeak(peaks_gr, TxDb=txdb,
                                      annoDb="org.Hs.eg.db")
            
            # Write annotation results
            annot_file <- file.path(output_dir, "tables", 
                                   sprintf("%s_peak_annotation.txt", condition))
            write.table(as.data.frame(peak_annot), annot_file,
                        sep="\t", quote=FALSE, row.names=FALSE)
            
            # Generate annotation summary
            annot_summary <- plotAnnoPie(peak_annot)
            pdf(file.path(output_dir, "plots", 
                         sprintf("%s_annotation_summary.pdf", condition)))
            print(annot_summary)
            dev.off()
            
            # Calculate overlap with V5 peaks
            log_message("Calculating peak overlaps", "INFO")
            overlaps <- findOverlaps(peaks_gr, v5_broad_peaks)
            n_overlaps <- length(unique(queryHits(overlaps)))
            
            # Write overlap statistics to report
            cat(sprintf("\nPeak Statistics for %s:\n", condition),
                sprintf("Total peaks: %d\n", length(peaks_gr)),
                sprintf("Peaks overlapping with V5: %d (%.1f%%)\n",
                        n_overlaps, 100 * n_overlaps/length(peaks_gr)),
                file=report_connection)
            
            # Initial H2AK119Ub analysis
            h2a_initial_stats <- generate_peak_stats(broad_peaks[[condition]], 
                                                   sprintf("Initial H2AK119Ub %s", condition))
            cat(h2a_initial_stats$text, file = report_connection)
            
            # Filter H2AK119Ub peaks
            # h2a_filtered <- filter_blacklist_regions(broad_peaks[[condition]], blacklist)
            h2a_filtered <- broad_peaks[[condition]]
            h2a_filtering_summary <- generate_filtering_summary(broad_peaks[[condition]], 
                                                             h2a_filtered, 
                                                             sprintf("H2AK119Ub %s Blacklist", condition))
            cat(h2a_filtering_summary, file = report_connection)
            
            # Generate and save overlap plots
            overlap_plots <- plot_overlap_analysis(v5_filtered, h2a_filtered, condition)
            
            if (!is.null(overlap_plots$width_dist)) {
                ggsave(file.path(plots_dir, 
                               sprintf("overlap_width_distribution_%s.pdf", 
                                     tolower(condition))),
                       overlap_plots$width_dist,
                       width = 8, height = 6)
            }
            
            ggsave(file.path(plots_dir, 
                           sprintf("overlap_summary_%s.pdf", 
                                 tolower(condition))),
                   overlap_plots$summary,
                   width = 8, height = 6)
            
            # Generate and save peak profile plot
            profile_plot <- plot_peak_profiles(v5_filtered, h2a_filtered, condition)
            if (!is.null(profile_plot)) {
                ggsave(file.path(plots_dir, 
                               sprintf("peak_profiles_%s.pdf", 
                                     tolower(condition))),
                       profile_plot,
                       width = 8, height = 6)
            }
            
            # Create combined plot for the condition
            plots_to_combine <- list(
                overlap_plots$width_dist,
                overlap_plots$summary
            )
            plots_to_combine <- plots_to_combine[!sapply(plots_to_combine, is.null)]
            
            if (length(plots_to_combine) > 0) {
                combined_plot <- do.call(grid.arrange, 
                                       c(plots_to_combine, 
                                         ncol = 1,
                                         top = sprintf("V5 vs H2AK119Ub %s Analysis", condition)))
                
                ggsave(file.path(plots_dir, 
                               sprintf("combined_analysis_%s.pdf", 
                                     tolower(condition))),
                       combined_plot,
                       width = 10, height = 12)
            }
        }
        
        log_message("Analysis completed successfully", "INFO", log_connection)
        cat("\nAnalysis completed successfully.\n", file = report_connection)
        
    }, error = function(e) {
        log_message(sprintf("ERROR: %s", e$message), "ERROR", log_connection)
        cat(sprintf("\nERROR: %s\n", e$message), file = report_connection)
        stop(e)
    })
}

# Set memory and performance optimization parameters
options(future.globals.maxSize = 2000 * 1024^2)  # 2GB limit for future package
options(mc.cores = parallel::detectCores() - 1)   # Use all but one core
options(scipen = 999)                           # Avoid scientific notation
options(datatable.optimize = 1)                 # Enable data.table optimization

# Function to plot peak profiles
plot_peak_profiles <- function(peaks, signal_data, condition, window_size = 2000) {
    # Calculate profile matrix
    profile_matrix <- do.call(cbind, lapply(signal_data, function(x) x$profile))
    
    # Calculate mean profile
    mean_profile <- rowMeans(profile_matrix, na.rm = TRUE)
    
    # Calculate standard error
    se_profile <- apply(profile_matrix, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    
    # Create position vector for x-axis
    positions <- seq(-window_size, window_size, length.out = length(mean_profile))
    
    # Create data frame for plotting
    plot_data <- data.frame(
        Position = positions,
        Mean = mean_profile,
        SE_plus = mean_profile + se_profile,
        SE_minus = mean_profile - se_profile
    )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = Position)) +
        geom_ribbon(aes(ymin = SE_minus, ymax = SE_plus), fill = "grey70", alpha = 0.3) +
        geom_line(aes(y = Mean), color = "blue", size = 1) +
        labs(x = "Distance from peak center (bp)",
             y = "Average signal",
             title = paste0(condition, " Peak Profiles")) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    
    # Save plot
    ggsave(file.path(output_dir, "plots", paste0(condition, "_peak_profiles.pdf")),
           p, width = 8, height = 6)
    
    return(p)
}

# Function to plot peak overlap analysis
plot_overlap_analysis <- function(v5_peaks, h2a_peaks, condition) {
    # Calculate overlaps
    overlaps <- findOverlaps(v5_peaks, h2a_peaks)
    
    # Calculate overlap statistics
    n_v5_peaks <- length(v5_peaks)
    n_h2a_peaks <- length(h2a_peaks)
    n_overlaps <- length(unique(queryHits(overlaps)))
    
    # Create Venn diagram data
    v5_name <- "V5 Peaks"
    h2a_name <- paste0(condition, " H2AK119Ub Peaks")
    
    # Set up colors
    colors <- brewer.pal(3, "Set2")[1:2]
    
    # Create Venn diagram
    venn_plot <- draw.pairwise.venn(
        area1 = n_v5_peaks,
        area2 = n_h2a_peaks,
        cross.area = n_overlaps,
        category = c(v5_name, h2a_name),
        lty = "blank",
        fill = colors,
        alpha = 0.5,
        cat.pos = c(0, 0),
        cat.dist = c(0.025, 0.025),
        cat.cex = 1,
        cex = 1.5,
        cat.col = colors,
        scaled = TRUE
    )
    
    # Save the plot
    pdf(file.path(output_dir, "plots", paste0(condition, "_overlap_analysis.pdf")),
        width = 8, height = 6)
    grid.draw(venn_plot)
    dev.off()
    
    # Return overlap statistics
    return(list(
        stats = list(
            v5_peaks = n_v5_peaks,
            h2a_peaks = n_h2a_peaks,
            overlapping_peaks = n_overlaps,
            percent_v5_overlap = 100 * n_overlaps / n_v5_peaks,
            percent_h2a_overlap = 100 * n_overlaps / n_h2a_peaks
        ),
        overlaps = overlaps
    ))
}

# Main execution block
main <- function() {
    # Enable garbage collection to manage memory
    gc()
    
    # Read V5 peaks with detailed tracking
    log_message("Starting to read V5 broad peaks...", "INFO")
    
    v5_broad_peaks_list <- list()
    for (rep_name in names(v5_broad_peaks_files)) {
        file_path <- v5_broad_peaks_files[[rep_name]]
        log_message(sprintf("Processing V5 replicate %s: %s", rep_name, file_path), "DEBUG")
        
        # Read peaks for this replicate
        peaks <- read_peaks(file_path, NULL)
        log_message(sprintf("Read %d peaks from replicate %s", length(peaks), rep_name), "DEBUG")
        
        v5_broad_peaks_list[[rep_name]] <- peaks
    }
    
    # Filter out any empty GRanges and merge V5 broad peaks from replicates
    valid_replicates <- v5_broad_peaks_list[sapply(v5_broad_peaks_list, length) > 0]
    log_message(sprintf("Number of valid V5 replicates: %d", length(valid_replicates)), "DEBUG")
    
    # Merge replicates
    v5_broad_peaks <- unlist(GRangesList(valid_replicates))
    log_message(sprintf("Total V5 peaks after merging replicates: %d", length(v5_broad_peaks)), "DEBUG")
    
    # Merge overlapping peaks with less aggressive parameters
    v5_broad_peaks <- merge_broad_peaks(v5_broad_peaks, 
                                      min_gap = 500,    # Reduced from 1000
                                      min_length = 500) # Reduced from 1000
    log_message(sprintf("Final number of V5 peaks after merging overlaps: %d", length(v5_broad_peaks)), "DEBUG")
    
    # Read H2AK119Ub broad peaks
    log_message("Reading H2AK119Ub broad peaks...", "INFO")
    broad_peaks <- read_peaks(h2a_peak_files, "broad")
    
    # Execute main analysis
    main_analysis(NULL, broad_peaks, NULL, v5_broad_peaks)
}

# Execute the main function
main()
