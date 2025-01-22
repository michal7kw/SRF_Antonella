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
})

# Set paths
h2a_base_dir <- "../SRF_H2AK119Ub/1_iterative_processing/analysis"
v5_peaks_file <- "../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_peaks.narrowPeak"
output_dir <- "results"

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# Function to standardize chromosome names
standardize_chromosomes <- function(gr) {
    # Get current seqlevels
    current_levels <- seqlevels(gr)
    
    # Add 'chr' prefix if missing and handle special cases
    new_levels <- current_levels
    for(i in seq_along(current_levels)) {
        if(!grepl("^chr", current_levels[i])) {
            if(current_levels[i] == "MT") {
                new_levels[i] <- "chrM"
            } else if(current_levels[i] == "23") {
                new_levels[i] <- "chrX"
            } else if(current_levels[i] == "24") {
                new_levels[i] <- "chrY"
            } else {
                new_levels[i] <- paste0("chr", current_levels[i])
            }
        }
    }
    
    # Create mapping
    names(new_levels) <- current_levels
    
    # Rename the seqlevels
    seqlevels(gr) <- new_levels
    
    # Keep only main chromosomes that exist in the data
    main_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
    existing_chroms <- intersect(seqlevels(gr), main_chroms)
    
    if(length(existing_chroms) == 0) {
        stop("No standard chromosomes found in the data")
    }
    
    gr <- keepSeqlevels(gr, existing_chroms, pruning.mode="coarse")
    return(gr)
}

# Function to analyze peak overlaps
analyze_overlaps <- function(peak_type) {
    message(paste("\nAnalyzing", peak_type, "peaks..."))
    
    # Read differential binding results
    diffbind_file <- file.path(h2a_base_dir, paste0("diffbind_", peak_type), "differential_peaks.csv")
    diff_peaks <- read.csv(diffbind_file)
    
    # Convert to GRanges
    diff_gr <- GRanges(
        seqnames = diff_peaks$seqnames,
        ranges = IRanges(start = diff_peaks$start, end = diff_peaks$end),
        strand = "*",
        score = diff_peaks$Fold,
        FDR = diff_peaks$FDR,
        Conc = diff_peaks$Conc
    )
    
    # Standardize chromosomes
    diff_gr <- standardize_chromosomes(diff_gr)
    
    # Read and standardize V5 peaks if not already done
    if(!exists("v5_peaks")) {
        message("Processing V5 peaks...")
        v5_peaks <<- import(v5_peaks_file)
        v5_peaks <<- standardize_chromosomes(v5_peaks)
    }
    
    # Ensure common chromosomes
    common_chroms <- intersect(seqlevels(diff_gr), seqlevels(v5_peaks))
    diff_gr <- keepSeqlevels(diff_gr, common_chroms, pruning.mode="coarse")
    v5_peaks_subset <- keepSeqlevels(v5_peaks, common_chroms, pruning.mode="coarse")
    
    # Find overlaps
    overlaps <- findOverlaps(diff_gr, v5_peaks_subset)
    
    # Create overlap statistics
    overlap_stats <- data.frame(
        Peak_Type = peak_type,
        Total_H2AK119Ub_Peaks = length(diff_gr),
        Total_V5_Peaks = length(v5_peaks_subset),
        Overlapping_Peaks = length(unique(queryHits(overlaps))),
        Overlap_Percentage = round(length(unique(queryHits(overlaps))) / length(diff_gr) * 100, 2)
    )
    
    # Analyze overlaps for significant peaks
    sig_peaks <- diff_gr[diff_gr$FDR < 0.05]
    sig_overlaps <- findOverlaps(sig_peaks, v5_peaks_subset)
    
    sig_stats <- data.frame(
        Peak_Type = paste(peak_type, "significant"),
        Total_H2AK119Ub_Peaks = length(sig_peaks),
        Total_V5_Peaks = length(v5_peaks_subset),
        Overlapping_Peaks = length(unique(queryHits(sig_overlaps))),
        Overlap_Percentage = round(length(unique(queryHits(sig_overlaps))) / length(sig_peaks) * 100, 2)
    )
    
    # Combine statistics
    stats <- rbind(overlap_stats, sig_stats)
    
    # Create Venn diagram
    venn_colors <- brewer.pal(3, "Set1")
    venn.plot <- venn.diagram(
        x = list(
            H2AK119Ub = 1:length(diff_gr),
            V5 = (length(diff_gr) + 1):(length(diff_gr) + length(v5_peaks_subset))
        ),
        category.names = c(paste0("H2AK119Ub\n(", peak_type, ")"), "V5"),
        filename = file.path(output_dir, "plots", paste0("venn_diagram_", peak_type, ".png")),
        output = TRUE,
        col = venn_colors[1:2],
        fill = alpha(venn_colors[1:2], 0.5),
        main = paste("Peak Overlap -", peak_type)
    )
    
    # Return statistics and overlapping peaks
    return(list(
        stats = stats,
        overlapping_peaks = diff_gr[queryHits(overlaps)],
        significant_overlapping_peaks = sig_peaks[queryHits(sig_overlaps)]
    ))
}

# Analyze both narrow and broad peaks
narrow_results <- analyze_overlaps("narrow")
broad_results <- analyze_overlaps("broad")

# Combine statistics
all_stats <- rbind(narrow_results$stats, broad_results$stats)
write.csv(all_stats, 
          file.path(output_dir, "tables", "overlap_statistics.csv"), 
          row.names = FALSE)

# Create comparison plot
ggplot(all_stats, aes(x = Peak_Type, y = Overlap_Percentage, fill = Peak_Type)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(title = "H2AK119Ub and V5 Peak Overlap Comparison",
         y = "Overlap Percentage",
         x = "Peak Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set1")
ggsave(file.path(output_dir, "plots", "overlap_comparison.pdf"), width = 8, height = 6)

# Analyze overlapping genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Function to get gene annotations
get_gene_annotations <- function(peaks, type) {
    peak_anno <- annotatePeak(peaks, 
                            TxDb = txdb,
                            tssRegion = c(-3000, 3000),
                            verbose = FALSE)
    
    genes_df <- as.data.frame(peak_anno)
    genes_df$peak_type <- type
    return(genes_df)
}

# Get annotations for all overlapping peaks
narrow_genes <- get_gene_annotations(narrow_results$overlapping_peaks, "narrow")
broad_genes <- get_gene_annotations(broad_results$overlapping_peaks, "broad")

# Combine gene annotations
all_genes <- rbind(narrow_genes, broad_genes)
write.csv(all_genes, 
          file.path(output_dir, "tables", "overlapping_genes.csv"), 
          row.names = FALSE)

# Create gene overlap Venn diagram
narrow_gene_symbols <- unique(narrow_genes$geneId)
broad_gene_symbols <- unique(broad_genes$geneId)

venn.diagram(
    x = list(Narrow = narrow_gene_symbols, Broad = broad_gene_symbols),
    filename = file.path(output_dir, "plots", "gene_overlap_venn.png"),
    category.names = c("Narrow Peaks", "Broad Peaks"),
    col = brewer.pal(3, "Set1")[1:2],
    fill = alpha(brewer.pal(3, "Set1")[1:2], 0.5),
    main = "Overlapping Genes Between Narrow and Broad Peaks"
)

# Create summary of results
summary_text <- c(
    "Cross-reference Analysis Summary",
    "==============================",
    "",
    "Peak Overlap Statistics:",
    paste("- Narrow peaks total:", narrow_results$stats$Total_H2AK119Ub_Peaks[1]),
    paste("- Narrow peaks overlapping with V5:", narrow_results$stats$Overlapping_Peaks[1], 
          sprintf("(%.1f%%)", narrow_results$stats$Overlap_Percentage[1])),
    paste("- Broad peaks total:", broad_results$stats$Total_H2AK119Ub_Peaks[1]),
    paste("- Broad peaks overlapping with V5:", broad_results$stats$Overlapping_Peaks[1],
          sprintf("(%.1f%%)", broad_results$stats$Overlap_Percentage[1])),
    "",
    "Gene Statistics:",
    paste("- Unique genes from narrow peaks:", length(narrow_gene_symbols)),
    paste("- Unique genes from broad peaks:", length(broad_gene_symbols)),
    paste("- Common genes between narrow and broad:", 
          length(intersect(narrow_gene_symbols, broad_gene_symbols)))
)

writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))

message("\nAnalysis completed successfully. Results are saved in the 'results' directory.")
