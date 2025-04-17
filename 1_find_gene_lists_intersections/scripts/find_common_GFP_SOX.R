#!/usr/bin/env Rscript

# Script to identify genes that have at least one peak of both GFP and SOX in their promoter regions
# Processes both standard and strict GFP peak calling results and compares them.

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(readr)
  library(GenomicFeatures)
})

# Function to read broadPeak files
read_broadpeak <- function(file_path) {
  cat("Parsing", file_path, "\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Error: File ", file_path, " does not exist")
  }
  
  # Read broadPeak file
  peaks <- read.table(file_path, header = FALSE, 
                     col.names = c("chrom", "start", "end", "name", "score", 
                                  "strand", "signalValue", "pValue", "qValue"))
  
  # Convert to GRanges object
  gr <- GRanges(
    seqnames = peaks$chrom,
    ranges = IRanges(start = peaks$start + 1, end = peaks$end),  # Convert to 1-based
    strand = "*",  # Peaks typically don't have strand information
    score = peaks$score,
    name = peaks$name
  )
  
  return(gr)
}

# Function to extract gene information from GTF file
extract_genes_from_gtf <- function(gtf_file) {
  cat("Extracting gene information from GTF file:", gtf_file, "\n")
  
  # Check if file exists
  if (!file.exists(gtf_file)) {
    stop("Error: File ", gtf_file, " does not exist")
  }
  
  # Import only gene features from GTF file to save memory
  genes <- import(gtf_file, feature.type="gene")
  
  # Extract essential information
  gene_info <- GRanges(
    seqnames = seqnames(genes),
    ranges = ranges(genes),
    strand = strand(genes),
    gene_id = genes$gene_id,
    gene_name = genes$gene_name
  )
  
  return(gene_info)
}

# Function to get promoter regions
get_promoters <- function(genes, promoter_size = 2000) {
  cat("Creating promoter regions with size:", promoter_size, "bp\n")
  
  # Create promoter regions based on strand
  promoters <- promoters(genes, upstream = promoter_size, downstream = 0)
  
  return(promoters)
}

# Function to find genes with peaks in promoter regions
find_genes_with_peaks_in_promoter <- function(peaks, promoters) {
  cat("Finding overlaps between peaks and promoters...\n")
  # Find overlaps between peaks and promoters
  overlaps <- findOverlaps(peaks, promoters)
  
  # Get unique gene IDs that have peaks in their promoters
  genes_with_peaks <- unique(promoters$gene_id[subjectHits(overlaps)])
  
  return(genes_with_peaks)
}

# Process data and save results function
process_data_and_save_results <- function(GFP_peaks, SOX_peaks, promoters, output_file) {
  # Find genes with peaks in promoter regions
  cat("Finding genes with GFP peaks in promoter regions\n")
  GFP_genes <- find_genes_with_peaks_in_promoter(GFP_peaks, promoters)
  cat("Found", length(GFP_genes), "genes with GFP peaks in promoter regions\n")
  
  cat("Finding genes with SOX peaks in promoter regions\n")
  SOX_genes <- find_genes_with_peaks_in_promoter(SOX_peaks, promoters)
  cat("Found", length(SOX_genes), "genes with SOX peaks in promoter regions\n")
  
  # Find genes with both GFP and SOX peaks in promoter regions
  common_genes <- intersect(GFP_genes, SOX_genes)
  cat("Found", length(common_genes), "genes with both GFP and SOX peaks in promoter regions\n")
  
  # Get gene information for common genes
  common_gene_info <- data.frame(
    gene_id = common_genes,
    gene_name = mcols(promoters[match(common_genes, promoters$gene_id)])$gene_name
  )
  
  # Sort by gene name
  common_gene_info <- common_gene_info %>% arrange(gene_name)
  
  # Save results to CSV
  write_csv(common_gene_info, output_file)
  cat("Results saved to", output_file, "\n")
  
  # Also save a BED file with the promoter regions of common genes for visualization
  bed_output_file <- sub("\\.csv$", "_promoters.bed", output_file)
  common_promoters <- promoters[promoters$gene_id %in% common_genes]
  
  # Convert to data frame for BED format
  bed_data <- data.frame(
    chrom = seqnames(common_promoters),
    start = start(common_promoters) - 1,  # Convert to 0-based for BED
    end = end(common_promoters),
    name = paste0(common_promoters$gene_name, "_", common_promoters$gene_id),
    score = 1000,
    strand = strand(common_promoters)
  )
  
  # Write BED file
  write.table(bed_data, bed_output_file, sep="\t", quote=FALSE, 
              row.names=FALSE, col.names=FALSE)
  cat("Promoter regions saved to", bed_output_file, "for visualization\n")
  
  return(common_gene_info)
}

# Main function
main <- function() {
  # Define base output directory
  output_dir_base <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/1_find_gene_lists_intersections/output"

  # File paths for inputs
  GFP_peak_file_standard <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling/GFP.broadPeak"
  GFP_peak_file_strict <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling_strict/GFP.broadPeak"
  SOX_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi_corrected/SOX2.broadPeak" # Ensure this path is correct for SOX2
  gtf_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"

  # File paths for outputs
  output_file_standard <- file.path(output_dir_base, "GFP_SOX_standard.csv")
  output_file_strict <- file.path(output_dir_base, "GFP_SOX_strict.csv")
  comparison_output <- file.path(output_dir_base, "GFP_SOX_comparison.csv")
  
  # Ensure output directory exists
  if (!dir.exists(output_dir_base)) {
    dir.create(output_dir_base, recursive = TRUE)
    cat("Created output directory:", output_dir_base, "\n")
  }

  # Extract gene information from GTF and get promoter regions (done once for all analyses)
  genes <- extract_genes_from_gtf(gtf_file)
  promoters <- get_promoters(genes)
  
  # Parse SOX peak file (done once for all analyses)
  SOX_peaks <- read_broadpeak(SOX_peak_file)
  
  # Process standard peaks
  cat("\n=== Processing Standard GFP vs SOX Peak Calling ===\n")
  GFP_peaks_standard <- read_broadpeak(GFP_peak_file_standard)
  standard_results <- process_data_and_save_results(GFP_peaks_standard, SOX_peaks, promoters, output_file_standard)
  
  # Process strict peaks
  cat("\n=== Processing Strict GFP vs SOX Peak Calling ===\n")
  GFP_peaks_strict <- read_broadpeak(GFP_peak_file_strict)
  strict_results <- process_data_and_save_results(GFP_peaks_strict, SOX_peaks, promoters, output_file_strict)
  
  # Find common genes between standard and strict GFP vs SOX analyses
  common_genes_comparison <- intersect(standard_results$gene_id, strict_results$gene_id)
  cat("\nFound", length(common_genes_comparison), "genes common to both standard and strict GFP vs SOX analyses\n")
  
  # Create a comparison file
  comparison_df <- data.frame(
    gene_id = common_genes_comparison,
    gene_name = standard_results$gene_name[match(common_genes_comparison, standard_results$gene_id)]
  ) %>% arrange(gene_name)
  
  write_csv(comparison_df, comparison_output)
  cat("Comparison results saved to", comparison_output, "\n")
}

# Run the main function
main()

# Remove the old TxDb implementation as it's not used in this multi-analysis structure
# main_txdb <- function() { ... }
