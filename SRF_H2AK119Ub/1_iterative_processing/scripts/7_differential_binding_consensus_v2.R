#' This script performs differential binding analysis on H2AK119Ub CUT&Tag data
#' comparing YAF and GFP samples, using pre-defined consensus peak sets for each condition.
#' It identifies regions with significantly different H2AK119Ub levels between conditions.
#'
#' Input files:
#' - analysis/6_consensus_peaks/GFP_consensus_peaks.bed: Consensus peaks for GFP condition.
#' - analysis/6_consensus_peaks/YAF_consensus_peaks.bed: Consensus peaks for YAF condition.
#' - analysis/3_alignment/{sample}.dedup.bam: Deduplicated alignment files for each sample.
#'
#' Output files:
#' - analysis/7_differential_binding_consensus_v2/: Directory containing DiffBind results.
#'
#' Dependencies:
#' - DiffBind, tidyverse, rtracklayer, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   org.Hs.eg.db, GenomicRanges, clusterProfiler, GenomeInfoDb

# Load required libraries
suppressPackageStartupMessages({
    library(DiffBind)
    library(tidyverse)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene) # For annotation if re-enabled
    library(org.Hs.eg.db) # For annotation if re-enabled
    library(GenomicRanges)
    library(clusterProfiler) # For annotation if re-enabled
    library(GenomeInfoDb)  # For renaming seqlevels
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("Usage: Rscript script.R OUTPUT_DIR ALIGN_DIR SAMPLE_IDS SAMPLE_CONDITIONS SAMPLE_REPLICATES GFP_CONSENSUS_PEAK_FILE YAF_CONSENSUS_PEAK_FILE")
}

# Define global constants from command line arguments
OUTPUT_DIR <- args[1]
ALIGN_DIR <- args[2]
# SAMPLE_IDS, SAMPLE_CONDITIONS, SAMPLE_REPLICATES are for BAM files
sample_ids_str <- args[3]
sample_conditions_str <- args[4]
sample_replicates_str <- args[5]
GFP_CONSENSUS_PEAK_FILE <- args[6]
YAF_CONSENSUS_PEAK_FILE <- args[7]

# ANNOTATION_DIR <- file.path(OUTPUT_DIR, "./annotation") # If annotation is re-enabled

# Parse sample information for BAM files
sample_ids <- strsplit(sample_ids_str, ",")[[1]]
sample_conditions <- strsplit(sample_conditions_str, ",")[[1]]
sample_replicates <- as.numeric(strsplit(sample_replicates_str, ",")[[1]])

if (length(sample_ids) != length(sample_conditions) || length(sample_ids) != length(sample_replicates)) {
    stop("Sample IDs, conditions, and replicates for BAM files must have the same length")
}

# Log a message with timestamp and level
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Check if all required input files exist
validate_files <- function(files_to_check) {
    missing_files <- files_to_check[!file.exists(files_to_check)]
    if (length(missing_files) > 0) {
        log_message(paste("Missing files:", paste(missing_files, collapse="\n")), level = "ERROR")
        stop("Missing files. Check the log for details.")
    }
    TRUE
}

# Create directories if they don't exist
create_dirs <- function(dirs_to_create) {
    for (d in dirs_to_create) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
}

# Standardize chromosome names in BED files to non-UCSC style (e.g., "1", "X", "MT")
# and ensure basic BED format compliance.
standardize_consensus_peak_file <- function(bed_file_path) {
    log_message(sprintf("Standardizing consensus peak file: %s", bed_file_path))
    tryCatch({
        if (!file.exists(bed_file_path) || file.info(bed_file_path)$size == 0) {
            stop(sprintf("File %s is missing or empty.", bed_file_path))
        }
        
        # Import BED file
        peaks_gr <- rtracklayer::import.bed(bed_file_path)
        
        if (length(peaks_gr) == 0) {
            log_message(sprintf("No peaks found in %s after import.bed. Check file format.", bed_file_path), level = "WARNING")
            # Write an empty valid BED file to avoid downstream errors if DiffBind requires a file
            # However, this scenario should ideally be caught by the shell script validation
            # For now, we'll let it proceed and DiffBind might error out, which is informative.
            # Or, write a minimal valid BED to a temp file and use that.
            # For now, we assume the shell script validation ensures non-empty files.
            # If it still results in zero peaks, DiffBind will handle it.
            return(bed_file_path) # Return original path, DiffBind will likely fail if empty
        }
        
        original_num_peaks <- length(peaks_gr)
        log_message(sprintf("Read %d peaks from %s", original_num_peaks, bed_file_path))

        # Standardize chromosome names: remove "chr" prefix
        seqlevels_style_orig <- seqlevelsStyle(peaks_gr)
        seqlevelsStyle(peaks_gr) <- "NCBI" # Converts "chr1" to "1", "chrX" to "X" etc.
                                         # If already "1", "X", it should remain as is.
        log_message(sprintf("Standardized seqlevels to NCBI style for %s. Original style: %s", bed_file_path, paste(seqlevels_style_orig, collapse=", ")))

        # Filter for standard chromosomes (non-UCSC format)
        standard_chroms_ncbi <- c(as.character(1:22), "X", "Y", "MT")
        original_seqlevels <- seqlevels(peaks_gr)
        peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "tidy", species = "Homo sapiens")
        
        current_num_peaks <- length(peaks_gr)
        log_message(sprintf("Filtered from %d to %d peaks based on standard NCBI chromosomes in %s", original_num_peaks, current_num_peaks, bed_file_path))
        
        if (current_num_peaks == 0) {
            stop(sprintf("No peaks left after standardizing and filtering chromosomes in %s. Original seqlevels: %s", bed_file_path, paste(original_seqlevels, collapse=", ")))
        }

        # Ensure essential BED columns are present (implicitly handled by GRanges to BED export)
        # rtracklayer::export.bed will create a valid BED file.
        # We overwrite the original file with the standardized one.
        rtracklayer::export.bed(peaks_gr, bed_file_path)
        log_message(sprintf("Successfully standardized and saved %d peaks to %s", current_num_peaks, bed_file_path))
        
    }, error = function(e) {
        log_message(sprintf("ERROR processing consensus peak file %s: %s", bed_file_path, e$message), level = "ERROR")
        stop(sprintf("Failed to process consensus peak file %s. Please check its format and content.", bed_file_path))
    })
    return(bed_file_path) # Return path to (potentially overwritten) standardized file
}


# Check and create BAM index files if missing
check_bam_indexes <- function(bam_files_paths) {
    for(bam_path in bam_files_paths) {
        index_path <- paste0(bam_path, ".bai")
        if(!file.exists(index_path)) {
            log_message(sprintf("BAM index %s not found. Creating index for %s", index_path, bam_path), "WARNING")
            # Ensure Rsamtools is available if using this. For now, rely on system samtools.
            # Rsamtools::indexBam(bam_path)
            system2("samtools", args = c("index", bam_path))
            if(!file.exists(index_path)) {
                stop(sprintf("Failed to create BAM index for %s", bam_path))
            }
            log_message(sprintf("Successfully created BAM index for %s", bam_path))
        }
    }
}

# --- Main Workflow ---

log_message("Starting differential binding analysis using consensus peaks.")

# Create output directories
# create_dirs(c(OUTPUT_DIR, ANNOTATION_DIR)) # ANNOTATION_DIR if re-enabled
create_dirs(c(OUTPUT_DIR))

# Define paths to BAM files
bam_read_paths <- file.path(ALIGN_DIR, paste0(sample_ids, ".dedup.bam"))

# Validate input files
files_to_validate <- c(bam_read_paths, GFP_CONSENSUS_PEAK_FILE, YAF_CONSENSUS_PEAK_FILE)
validate_files(files_to_validate)
log_message("All specified input files exist.")

# Check/Create BAM indexes
check_bam_indexes(bam_read_paths)
log_message("BAM indexes checked/created.")

# Standardize consensus peak files
log_message("Standardizing consensus peak files...")
standardize_consensus_peak_file(GFP_CONSENSUS_PEAK_FILE)
standardize_consensus_peak_file(YAF_CONSENSUS_PEAK_FILE)
log_message("Consensus peak files standardized.")

# Create samples data frame for DiffBind
# Peaks column will point to the *same* consensus peak file for all samples of that condition.
peak_files_for_samples <- ifelse(sample_conditions == "GFP", GFP_CONSENSUS_PEAK_FILE, YAF_CONSENSUS_PEAK_FILE)

samples_df <- data.frame(
    SampleID = sample_ids,
    Factor = rep("H2AK119Ub", length(sample_ids)), # Assuming Factor is constant
    Condition = sample_conditions,
    Replicate = sample_replicates,
    bamReads = bam_read_paths,
    Peaks = peak_files_for_samples,
    PeakCaller = rep("bed", length(sample_ids)) # Using BED format for consensus peaks
)

log_message("Sample sheet for DiffBind:")
print(samples_df)

# Save the constructed sample sheet
write.csv(samples_df, file.path(OUTPUT_DIR, "consensus_sample_sheet.csv"), row.names = FALSE)

# Create DiffBind object
log_message("Creating DiffBind object using consensus peaks...")
# minOverlap = 1 means the peak must be in the provided (consensus) peak set.
# Since we provide a single consensus peak file per condition, this effectively means
# DiffBind will create a master peak set from the union of GFP_consensus and YAF_consensus peaks.
dba_data <- dba(sampleSheet = samples_df,
                minOverlap = 1,
                peakFormat = "bed", # Peak files are in BED format
                config = data.frame(AnalysisMethod = "DBA_DESEQ2",
                                    fragmentSize = 150, # Adjust if necessary
                                    doBlacklist = FALSE, # Assuming consensus peaks are already filtered
                                    RunParallel = TRUE,
                                    summits = FALSE)) # No summits for broad/consensus regions

log_message("Initial peak counts per sample (based on consensus sets):")
# This will show how many peaks from the *master set* (union of consensus) overlap with each sample's *own* consensus set.
# Since Peaks column points to the condition's consensus, it will show full count for that condition's consensus.
print(dba.show(dba_data))

# Count reads
log_message("Counting reads in consensus peak regions...")
tryCatch({
    dba_data <- dba.count(dba_data,
                          bUseSummarizeOverlaps = TRUE,
                          filter = 0, 
                          minCount = 0,
                          bRemoveDuplicates = TRUE, # From BAMs
                          score = DBA_SCORE_READS, # Count reads
                          summits = FALSE) # No summits for broad/consensus regions
    
    log_message("dba.count completed. After counting - sample statistics:")
    count_info <- dba.show(dba_data)
    print(count_info)
    
    if(is.null(count_info) || sum(count_info$Reads, na.rm = TRUE) == 0) {
         log_message("WARNING: No reads were counted across any samples, or count_info is NULL/zero.", level="WARNING")
    }
    
    log_message("Full dba_data object state after dba.count:")
    print(dba_data)
    
    if (is.null(dba_data$binding) || nrow(dba_data$binding) == 0) {
        log_message("CRITICAL: dba_data$binding (counts table) is NULL or has 0 rows after dba.count.", level="ERROR")
        stop("CRITICAL: Counts table (dba_data$binding) is empty after dba.count.")
    }

}, error = function(e) {
    log_message(sprintf("ERROR during dba.count: %s", e$message), level = "ERROR")
    log_message("Printing dba_data object state before dba.count error:", level = "ERROR")
    print(dba_data)
    stop(sprintf("Failed during dba.count: %s", e$message))
})

# Normalize and perform differential analysis
log_message("Performing differential analysis...")
dba_data <- dba.normalize(dba_data, 
                      normalize = DBA_NORM_LIB, # Or other appropriate normalization
                      background = FALSE) # No background normalization if no control samples

log_message("Setting up contrasts...")
# Ensure Condition is a factor for dba.contrast
dba_data$samples$Condition <- factor(dba_data$samples$Condition)
dba_data <- dba.contrast(dba_data, 
                         categories = DBA_CONDITION,
                         minMembers = 2) # Ensure at least 2 replicates per condition for comparison

log_message("Contrast information:")
print(dba.show(dba_data, bContrasts = TRUE))

tryCatch({
    log_message("Inspecting dba_data object just before dba.analyze call:")
    if(is.null(dba_data$binding) || nrow(dba_data$binding) == 0) {
        log_message("ERROR: dba_data$binding (counts) is null or has 0 rows before dba.analyze!", level="ERROR")
    } else {
        log_message(sprintf("Number of sites (rows in count table) in dba_data: %d", nrow(dba_data$binding)))
        log_message(sprintf("Number of samples (cols in count table) in dba_data: %d", ncol(dba_data$binding)))
    }
    
    log_message("Starting DESeq2 analysis with dispersion estimation...")
    dba_data <- dba.analyze(dba_data, method = DBA_DESEQ2)
    log_message("DESeq2 analysis completed successfully")
    
    analysis_stats <- dba.show(dba_data, bContrasts = TRUE)
    log_message("Analysis statistics:")
    print(analysis_stats)
    
    # Generate plots
    pdf(file.path(OUTPUT_DIR, "consensus_diffbind_plots.pdf"))
    tryCatch({
        dba.plotHeatmap(dba_data, contrast=1) # Assuming first contrast
        title("Sample correlation heatmap (Consensus Peaks)")
        
        dba.plotPCA(dba_data, contrast=1, attributes=DBA_CONDITION, label=DBA_ID)
        title("PCA plot of samples (Consensus Peaks)")
        
        dba.plotMA(dba_data, contrast=1)
        title("MA plot of differential binding (Consensus Peaks)")
        
        dba.plotVolcano(dba_data, contrast=1)
        title("Volcano plot of differential binding (Consensus Peaks)")
    }, error = function(e_plot) {
        log_message(sprintf("Warning: Could not generate some plots: %s", e_plot$message), "WARNING")
    }, finally = {
        dev.off()
    })
    
}, error = function(e_analyze) {
    log_message(sprintf("Error in differential analysis: %s", e_analyze$message), "ERROR")
    stop(sprintf("Differential analysis failed: %s", e_analyze$message))
})

log_message("Extracting results...")
# th=1 means all sites are reported, FDR/p-value can be used for filtering later
dba_results_gr <- dba.report(dba_data,
                         th = 1,      # Report all sites
                         bCalled = FALSE, # Not relevant for consensus approach like this
                         bNormalized = TRUE,
                         bCounts = TRUE,  # Include normalized and raw counts
                         contrast = 1) # Assuming the first contrast is the one of interest

if (is.null(dba_results_gr) || length(dba_results_gr) == 0) {
    log_message("No differential binding sites found or dba.report returned NULL.", "WARNING")
    df_results <- data.frame() # Create an empty data frame
} else {
    log_message(sprintf("Extracted %d differential binding sites.", length(dba_results_gr)))
    df_results <- as.data.frame(dba_results_gr)
}


log_message("Saving results...")
saveRDS(dba_data, file.path(OUTPUT_DIR, "dba_consensus_analysis.rds"))
saveRDS(dba_results_gr, file.path(OUTPUT_DIR, "dba_consensus_significant_peaks.rds"))
write.csv(df_results, file.path(OUTPUT_DIR, "dba_consensus_differential_peaks.csv"), row.names = FALSE)

log_message(sprintf("Found %d differential binding sites (all reported, th=1).", nrow(df_results)))

# Optional: Peak Annotation (can be re-enabled if needed)
# log_message("Peak annotation step is currently disabled.")
# if (nrow(df_results) > 0 & !is.null(dba_results_gr) & length(dba_results_gr) > 0) {
#     log_message("Annotating peaks...")
#     # For annotation, UCSC-style names might be needed by TxDb.Hsapiens.UCSC.hg38.knownGene
#     # The dba_results_gr should have NCBI style from standardization. Convert if needed.
#     peaks_to_annotate <- dba_results_gr
#     seqlevelsStyle(peaks_to_annotate) <- "UCSC" # Convert "1" to "chr1"
#
#     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#     peakAnno <- annotatePeak(peaks_to_annotate, TxDb = txdb,
#                              tssRegion = c(-3000, 3000),
#                              verbose = FALSE, annoDb="org.Hs.eg.db")
#
#     dir.create(ANNOTATION_DIR, showWarnings = FALSE, recursive = TRUE)
#     saveRDS(peakAnno, file.path(ANNOTATION_DIR, "consensus_peak_annotation.rds"))
#     write.csv(as.data.frame(peakAnno),
#               file.path(ANNOTATION_DIR, "consensus_peak_annotation.csv"),
#               row.names = FALSE)
#
#     pdf(file.path(ANNOTATION_DIR, "consensus_annotation_plots.pdf"))
#     plotAnnoPie(peakAnno)
#     plotDistToTSS(peakAnno)
#     dev.off()
#     log_message("Peak annotation completed.")
# } else {
#     log_message("Skipping peak annotation as no significant peaks were found or results were empty.")
# }

log_message("Successfully completed differential binding analysis with consensus peaks.")