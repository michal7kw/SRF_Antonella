# Load required libraries
suppressPackageStartupMessages({
    library(DiffBind)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(clusterProfiler)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
})

# Logging function
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Input validation
validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
    TRUE
}

# Create output directories
create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Standardize chromosome names
standardize_chromosomes <- function(gr) {
    # Add 'chr' prefix if missing
    current_chroms <- seqlevels(gr)
    if (!all(grepl("^chr", current_chroms))) {
        new_levels <- paste0("chr", current_chroms)
        names(new_levels) <- current_chroms
        seqlevels(gr) <- new_levels
    }
    
    # Handle mitochondrial chromosome naming
    current_chroms <- seqlevels(gr)
    mt_matches <- grep("^chr(M|MT)$", current_chroms, value = TRUE)
    if (length(mt_matches) > 0) {
        # Rename to standard chrM
        new_levels <- current_chroms
        names(new_levels) <- current_chroms
        new_levels[mt_matches] <- "chrM"
        seqlevels(gr) <- new_levels
    }
    
    # Keep only standard chromosomes
    standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    gr <- keepSeqlevels(gr, standard_chroms, pruning.mode = "coarse")
    
    # Set genome
    genome(gr) <- "hg38"
    
    return(gr)
}

# Save plots with error handling
save_plot <- function(plot, filename, width = 8, height = 6) {
    tryCatch({
        pdf(filename, width = width, height = height)
        print(plot)
        dev.off()
        log_message(sprintf("Saved plot to %s", filename))
    }, error = function(e) {
        log_message(sprintf("Error saving plot to %s: %s", filename, e$message), "ERROR")
    })
}

# Create MA plot
create_ma_plot <- function(data) {
    ggplot(data, aes(x = log2(Conc), y = Fold)) +
        geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
        scale_color_manual(values = c("grey50", "red3")) +
        scale_size_continuous(range = c(0.5, 2)) +
        labs(x = "log2 Mean Concentration",
             y = "log2 Fold Change",
             title = "MA Plot") +
        theme_minimal()
}

# Create volcano plot
create_volcano_plot <- function(data) {
    ggplot(data, aes(x = Fold, y = -log10(FDR))) +
        geom_point(aes(color = FDR < 0.05, size = abs(Fold)), alpha = 0.6) +
        scale_color_manual(values = c("grey50", "red3")) +
        scale_size_continuous(range = c(0.5, 2)) +
        labs(x = "log2 Fold Change",
             y = "-log10(FDR)",
             title = "Volcano Plot") +
        theme_minimal()
}

# Calculate summary statistics
calculate_summary_stats <- function(data) {
    data.frame(
        Metric = c("Total Peaks", 
                  "Significant Peaks", 
                  "YAF-enriched", 
                  "GFP-enriched"),
        Count = c(
            nrow(data),
            sum(data$FDR < 0.05),
            sum(data$FDR < 0.05 & data$Fold > 0),
            sum(data$FDR < 0.05 & data$Fold < 0)
        )
    )
} 