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
})

# Set paths
h2a_base_dir <- "../SRF_H2AK119Ub/1_iterative_processing/analysis"
v5_peaks_file <- "../SRF_V5/peaks/SES-V5ChIP-Seq2_S6_peaks.narrowPeak"
output_dir <- "results"

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

# Add distance-based analysis
analyze_peak_distances <- function(h2a_peaks, v5_peaks) {
    # Calculate distances between H2AK119Ub and nearest V5 peaks
    distances <- distanceToNearest(h2a_peaks, v5_peaks)
    
    # Create distance distribution plot
    dist_plot <- ggplot(data.frame(distance=mcols(distances)$distance), 
                       aes(x=distance)) +
        geom_histogram(binwidth=100) +
        labs(title="Distance Distribution between H2AK119Ub and V5 Peaks",
             x="Distance (bp)", y="Count")
    
    return(list(distances=distances, plot=dist_plot))
}

# Add signal correlation analysis
analyze_signal_correlation <- function(h2a_peaks, v5_peaks, h2a_bigwig, v5_bigwig) {
    # Import signal data from bigWig files
    h2a_signal <- import(h2a_bigwig, as="RleList")
    v5_signal <- import(v5_bigwig, as="RleList")
    
    # Extract signal values at peaks
    h2a_matrix <- matrix(nrow=length(h2a_peaks), ncol=2)
    v5_matrix <- matrix(nrow=length(v5_peaks), ncol=2)
    
    # Calculate correlation
    for(i in seq_along(h2a_peaks)) {
        h2a_matrix[i,] <- c(
            mean(h2a_signal[[as.character(seqnames(h2a_peaks[i]))]][start(h2a_peaks[i]):end(h2a_peaks[i])]),
            mean(v5_signal[[as.character(seqnames(h2a_peaks[i]))]][start(h2a_peaks[i]):end(h2a_peaks[i])])
        )
    }
    
    # Generate correlation plot
    correlation_plot <- ggplot(data.frame(H2A=h2a_matrix[,1], V5=h2a_matrix[,2]), 
                             aes(x=H2A, y=V5)) +
        geom_point(alpha=0.5) +
        geom_smooth(method="lm") +
        labs(title="Signal Correlation between H2AK119Ub and V5",
             x="H2AK119Ub Signal", y="V5 Signal") +
        theme_bw()
    
    # Calculate correlation coefficient
    cor_coef <- cor(h2a_matrix[,1], h2a_matrix[,2], method="pearson")
    
    return(list(correlation=cor_coef, plot=correlation_plot, signal_values=h2a_matrix))
}

# Add statistical significance testing
assess_overlap_significance <- function(peak1, peak2, genome_size, n_permutations=1000) {
    observed_overlap <- length(findOverlaps(peak1, peak2))
    
    # Perform permutation test
    null_distribution <- numeric(n_permutations)
    for(i in 1:n_permutations) {
        # Randomly shuffle peaks
        shuffled_peaks <- shuffle(peak2, genome=genome_size)
        null_distribution[i] <- length(findOverlaps(peak1, shuffled_peaks))
    }
    
    # Calculate empirical p-value
    p_value <- sum(null_distribution >= observed_overlap) / n_permutations
    
    return(list(p_value=p_value, null_distribution=null_distribution))
}

# Add genomic feature enrichment analysis
analyze_genomic_context <- function(peaks, txdb) {
    # Analyze enrichment in genomic features
    genomic_annotation <- annotatePeak(peaks,
                                     TxDb=txdb,
                                     tssRegion=c(-3000, 3000),
                                     verbose=FALSE)
    
    # Create feature distribution plot
    feature_plot <- plotAnnoPie(genomic_annotation)
    
    # Perform enrichment analysis
    feature_enrichment <- ChIPseeker::enrichPeakOverlap(peaks,
                                                       TxDb=txdb,
                                                       nShuffle=1000,
                                                       chainFile=NULL)
    
    # Create enrichment plot
    enrichment_plot <- ggplot(data.frame(feature_enrichment), 
                            aes(x=reorder(feature, log2FE), y=log2FE)) +
        geom_bar(stat="identity") +
        coord_flip() +
        theme_bw() +
        labs(title="Genomic Feature Enrichment",
             x="Feature", y="Log2 Fold Enrichment")
    
    return(list(
        annotation=genomic_annotation,
        enrichment=feature_enrichment,
        plots=list(feature_dist=feature_plot, 
                  enrichment=enrichment_plot)
    ))
}

# Add QC metrics
calculate_qc_metrics <- function(peaks, genome_size) {
    # Calculate peak width distribution
    widths <- width(peaks)
    width_stats <- list(
        mean=mean(widths),
        median=median(widths),
        sd=sd(widths),
        quantiles=quantile(widths, probs=c(0.25, 0.75))
    )
    
    # Calculate genome coverage
    coverage <- sum(widths) / genome_size
    
    # Calculate peak score distribution
    scores <- peaks$score
    score_stats <- list(
        mean=mean(scores),
        median=median(scores),
        sd=sd(scores),
        quantiles=quantile(scores, probs=c(0.25, 0.75))
    )
    
    # Create QC plots
    width_plot <- ggplot(data.frame(width=widths), aes(x=width)) +
        geom_histogram(bins=50) +
        theme_bw() +
        labs(title="Peak Width Distribution",
             x="Width (bp)", y="Count")
    
    score_plot <- ggplot(data.frame(score=scores), aes(x=score)) +
        geom_histogram(bins=50) +
        theme_bw() +
        labs(title="Peak Score Distribution",
             x="Score", y="Count")
    
    return(list(
        width_stats=width_stats,
        coverage=coverage,
        score_stats=score_stats,
        plots=list(width=width_plot, score=score_plot)
    ))
}

# Add pathway analysis
perform_pathway_analysis <- function(gene_list) {
    # Convert gene IDs if necessary
    if(!all(grepl("^ENTREZ:", gene_list))) {
        gene_list <- bitr(gene_list, 
                         fromType="SYMBOL",
                         toType="ENTREZID",
                         OrgDb=org.Hs.eg.db)$ENTREZID
    }
    
    # Perform GO enrichment
    go_enrichment <- enrichGO(gene=gene_list,
                             OrgDb=org.Hs.eg.db,
                             ont="ALL",
                             pAdjustMethod="BH",
                             pvalueCutoff=0.05,
                             qvalueCutoff=0.2)
    
    # Perform KEGG pathway analysis
    kegg_enrichment <- enrichKEGG(gene=gene_list,
                                 organism="hsa",
                                 pAdjustMethod="BH",
                                 pvalueCutoff=0.05,
                                 qvalueCutoff=0.2)
    
    # Create enrichment plots
    go_plot <- dotplot(go_enrichment, showCategory=20) +
        ggtitle("GO Term Enrichment")
    
    kegg_plot <- dotplot(kegg_enrichment, showCategory=20) +
        ggtitle("KEGG Pathway Enrichment")
    
    return(list(
        go=go_enrichment,
        kegg=kegg_enrichment,
        plots=list(go=go_plot, kegg=kegg_plot)
    ))
}

# Execute additional analyses
message("\nPerforming additional analyses...")

# Perform distance analysis
narrow_distances <- analyze_peak_distances(narrow_results$overlapping_peaks, v5_peaks)
broad_distances <- analyze_peak_distances(broad_results$overlapping_peaks, v5_peaks)

# Save distance plots
ggsave(file.path(output_dir, "plots", "narrow_peak_distances.pdf"), 
       narrow_distances$plot, width=8, height=6)
ggsave(file.path(output_dir, "plots", "broad_peak_distances.pdf"), 
       broad_distances$plot, width=8, height=6)

# Perform genomic context analysis
narrow_context <- analyze_genomic_context(narrow_results$overlapping_peaks, txdb)
broad_context <- analyze_genomic_context(broad_results$overlapping_peaks, txdb)

# Save genomic context plots
pdf(file.path(output_dir, "plots", "genomic_context.pdf"))
grid.arrange(narrow_context$plots$feature_dist, broad_context$plots$feature_dist,
            narrow_context$plots$enrichment, broad_context$plots$enrichment,
            ncol=2)
dev.off()

# Calculate QC metrics
narrow_qc <- calculate_qc_metrics(narrow_results$overlapping_peaks, 3.2e9)  # Human genome size
broad_qc <- calculate_qc_metrics(broad_results$overlapping_peaks, 3.2e9)

# Save QC plots
pdf(file.path(output_dir, "plots", "qc_metrics.pdf"))
grid.arrange(narrow_qc$plots$width, broad_qc$plots$width,
            narrow_qc$plots$score, broad_qc$plots$score,
            ncol=2)
dev.off()

# Perform pathway analysis for overlapping genes
narrow_pathways <- perform_pathway_analysis(narrow_gene_symbols)
broad_pathways <- perform_pathway_analysis(broad_gene_symbols)

# Save pathway plots
pdf(file.path(output_dir, "plots", "pathway_analysis.pdf"))
grid.arrange(narrow_pathways$plots$go, broad_pathways$plots$go,
            narrow_pathways$plots$kegg, broad_pathways$plots$kegg,
            ncol=2)
dev.off()

message("Additional analyses completed. Results saved in the output directory.")

# Perform signal correlation analysis
message("Performing signal correlation analysis...")

# Calculate average bigwig signal for H2AK119Ub
h2a_signal_correlation <- analyze_signal_correlation(
    narrow_results$overlapping_peaks,
    v5_peaks,
    h2a_bigwig_files$YAF$rep1,  # You might want to average all replicates
    v5_bigwig_files$ChIP
)

# Save correlation plot
ggsave(file.path(output_dir, "plots", "signal_correlation.pdf"),
       h2a_signal_correlation$plot,
       width = 8, height = 6)

# Add correlation results to summary
cat("\nSignal Correlation Analysis:",
    paste("Pearson correlation coefficient:", round(h2a_signal_correlation$correlation, 3)),
    file = file.path(output_dir, "analysis_summary.txt"),
    append = TRUE,
    sep = "\n")
