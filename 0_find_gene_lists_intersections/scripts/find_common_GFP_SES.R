#!/usr/bin/env Rscript

# Script to identify genes that have at least one peak of both GFP and SES in their promoter regions
# This is a more efficient R implementation using GenomicRanges

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(readr)
  library(GenomicFeatures)
})

# BASE_PATH <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
BASE_PATH <- "/mnt/d/Github/SRF_H2AK119Ub_cross_V5"

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
# Function to read BED files
read_bed <- function(file_path) {
  cat("Parsing BED file:", file_path, "\n")
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Error: File ", file_path, " does not exist")
  }
  
  # Import BED file using rtracklayer
  # This function automatically handles 0-based vs 1-based and creates a GRanges object
  gr <- rtracklayer::import(file_path, format = "bed")
  cat("Successfully parsed", file_path, "- found", length(gr), "regions\n")
  cat("Original seqnames sample from BED", basename(file_path), ":", paste(head(unique(as.character(seqnames(gr)))), collapse=", "), "\n")
  tryCatch({
    current_style <- seqlevelsStyle(gr)
    # Attempt to set to UCSC if not already a recognized style or if it's NCBI/Ensembl (Gencode GTF is UCSC)
    if (!("UCSC" %in% current_style)) {
        cat("Attempting to set seqlevels style to UCSC for", basename(file_path), "(current: ", paste(current_style, collapse="/"), ")\n")
        seqlevelsStyle(gr) <- "UCSC"
        cat("New seqnames sample after UCSC attempt:", paste(head(unique(as.character(seqnames(gr)))), collapse=", "), "\n")
    } else {
        cat("Seqlevels style for", basename(file_path), "is already UCSC. No change applied.\n")
    }
  }, warning = function(w) {
    cat("Warning while setting/checking seqlevels style to UCSC for", basename(file_path), ":", conditionMessage(w), "\n")
  }, error = function(e) {
    cat("Error while setting/checking seqlevels style to UCSC for", basename(file_path), ":", conditionMessage(e), "\n")
    cat("Proceeding with original seqlevels for", basename(file_path), "\n")
  })
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
  cat("Original seqnames sample from GTF:", paste(head(unique(as.character(seqnames(genes)))), collapse=", "), "\n")
  tryCatch({
    current_style <- seqlevelsStyle(genes)
    if (!("UCSC" %in% current_style)) { # Gencode is UCSC, this ensures it if auto-detection failed
        cat("Attempting to set seqlevels style to UCSC for GTF genes (current: ", paste(current_style, collapse="/"), ")\n")
        seqlevelsStyle(genes) <- "UCSC"
        cat("New seqnames sample for GTF after UCSC attempt:", paste(head(unique(as.character(seqnames(genes)))), collapse=", "), "\n")
    } else {
        cat("Seqlevels style for GTF is already UCSC. No change applied.\n")
    }
  }, warning = function(w) {
    cat("Warning while setting/checking seqlevels style to UCSC for GTF genes:", conditionMessage(w), "\n")
  }, error = function(e) {
    cat("Error while setting/checking seqlevels style to UCSC for GTF genes:", conditionMessage(e), "\n")
    cat("Proceeding with original seqlevels for GTF genes\n")
  })
  
  # Extract essential information
  gene_info <- GRanges( # This will inherit seqinfo from 'genes'
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
process_data_and_save_results <- function(GFP_peaks, SES_peaks, promoters, output_file) {
  # Find genes with peaks in promoter regions
  cat("Finding genes with GFP peaks in promoter regions\n")
  GFP_genes <- find_genes_with_peaks_in_promoter(GFP_peaks, promoters)
  cat("Found", length(GFP_genes), "genes with GFP peaks in promoter regions\n")
  
  cat("Finding genes with SES peaks in promoter regions\n")
  SES_genes <- find_genes_with_peaks_in_promoter(SES_peaks, promoters)
  cat("Found", length(SES_genes), "genes with SES peaks in promoter regions\n")
  
  # Find genes with both GFP and SES peaks in promoter regions
  common_genes <- intersect(GFP_genes, SES_genes)
  cat("Found", length(common_genes), "genes with both GFP and SES peaks in promoter regions\n")
  
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
  
  if (length(common_promoters) > 0) {
    # Convert to data frame for BED format
    bed_data <- data.frame(
      chrom = as.character(seqnames(common_promoters)), # Ensure character to avoid factor issues
      start = start(common_promoters) - 1,  # Convert to 0-based for BED
      end = end(common_promoters),
      name = paste0(common_promoters$gene_name, "_", common_promoters$gene_id),
      score = 1000, # Or some other relevant score
      strand = as.character(strand(common_promoters)) # Ensure character
    )
    
    # Write BED file
    write.table(bed_data, bed_output_file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    cat("Promoter regions saved to", bed_output_file, "for visualization\n")
  } else {
    cat("No common promoter regions to save to BED file for", basename(output_file), "(0 common genes found)\n")
    # To prevent downstream tools from failing on an empty file, one might write a header-only BED or an empty file.
    # For now, it just prints a message and doesn't create the BED.
    # writeLines(c("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand"), bed_output_file) # Example: header only
  }
  
  return(common_gene_info)
}

# Main function
main <- function() {
  # File paths
  GFP_peak_file_standard <- file.path(BASE_PATH, "SRF_H2AK119Ub", "1_iterative_processing", "analysis", "6_consensus_peaks", "GFP_consensus_peaks.bed")
  output_file_standard <- file.path(BASE_PATH, "0_find_gene_lists_intersections/output", "GFP_SES.csv")

  GFP_peak_file_strict <- file.path(BASE_PATH, "SRF_H2AK119Ub", "1_iterative_processing", "analysis", "6_consensus_peaks_strict", "GFP_consensus_peaks.bed")
  output_file_strict <- file.path(BASE_PATH, "0_find_gene_lists_intersections/output", "GFP_SES_strict.csv")
  
  SES_peak_file <- file.path(BASE_PATH, "SRF_SES_V5/results_data_from_ncbi_corrected/SES.broadPeak")
  gtf_file <- file.path(BASE_PATH, "COMMON_DATA/gencode.v43.basic.annotation.gtf")
  
  # Extract gene information from GTF and get promoter regions (done once for both analyses)
  genes <- extract_genes_from_gtf(gtf_file)
  promoters <- get_promoters(genes)
  
  # Parse SES peak file (done once for both analyses)
  SES_peaks <- read_broadpeak(SES_peak_file)
  
  # Process standard peaks
  cat("\n=== Processing Standard Peak Calling ===\n")
  GFP_peaks_standard <- read_bed(GFP_peak_file_standard) # Use read_bed for .bed file
  standard_results <- process_data_and_save_results(GFP_peaks_standard, SES_peaks, promoters, output_file_standard)
  
  # Process strict peaks
  cat("\n=== Processing Strict Peak Calling ===\n")
  GFP_peaks_strict <- read_bed(GFP_peak_file_strict) # Use read_bed for .bed file
  strict_results <- process_data_and_save_results(GFP_peaks_strict, SES_peaks, promoters, output_file_strict)
  
  # Find common genes between standard and strict analyses
  common_genes <- intersect(standard_results$gene_id, strict_results$gene_id)
  cat("\nFound", length(common_genes), "genes common to both standard and strict analyses\n")
  
  # Create a comparison file
  comparison_output <- file.path(BASE_PATH, "0_find_gene_lists_intersections/output", "GFP_SES_comparison.csv") # Corrected path
  comparison_df <- data.frame(
    gene_id = common_genes,
    gene_name = standard_results$gene_name[match(common_genes, standard_results$gene_id)]
  ) %>% arrange(gene_name)
  
  write_csv(comparison_df, comparison_output)
  cat("Comparison results saved to", comparison_output, "\n")
}

# Run the main function
main()