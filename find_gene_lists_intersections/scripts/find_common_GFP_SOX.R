#!/usr/bin/env Rscript

# Script to identify genes that have at least one peak of both GFP and SOX in their promoter regions
# This is a more efficient R implementation using GenomicRanges

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

# Main function
main <- function() {
  # File paths

  GFP_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling/GFP.broadPeak"
  output_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/GFP_SOX_strict.csv"

  # GFP_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling/GFP.broadPeak"
  # output_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/GFP_SOX.csv"
  
  SOX_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi_corrected/SOX2.broadPeak"
  gtf_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"
  
  # Parse peak files
  GFP_peaks <- read_broadpeak(GFP_peak_file)
  SOX_peaks <- read_broadpeak(SOX_peak_file)
  
  # Extract gene information from GTF
  genes <- extract_genes_from_gtf(gtf_file)
  
  # Get promoter regions
  promoters <- get_promoters(genes)
  
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
}

# Alternative implementation using TxDb package (faster for large GTF files)
main_txdb <- function() {
  cat("Using TxDb approach (alternative implementation)\n")
  
  # File paths

  GFP_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling/GFP.broadPeak"
  output_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/GFP_SOX_strict.csv"

  # GFP_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling/GFP.broadPeak"
  # output_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/GFP_SOX.csv"
  
  SOX_peak_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_SES_V5/results_data_from_ncbi_corrected/SOX2.broadPeak"
  gtf_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"
  
  # Parse peak files
  GFP_peaks <- read_broadpeak(GFP_peak_file)
  SOX_peaks <- read_broadpeak(SOX_peak_file)
  
  # Create TxDb object from GTF file
  cat("Creating TxDb from GTF file...\n")
  txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
  
  # Get gene information
  cat("Extracting gene information...\n")
  genes <- genes(txdb)
  
  # Add gene names from GTF (since TxDb doesn't store gene names)
  gtf_genes <- import(gtf_file, feature.type="gene")
  gene_names <- gtf_genes$gene_name
  names(gene_names) <- gtf_genes$gene_id
  
  # Create promoter regions
  cat("Creating promoter regions...\n")
  promoters <- promoters(genes, upstream=2000, downstream=0)
  
  # Find genes with peaks in promoter regions
  cat("Finding genes with GFP peaks in promoter regions\n")
  GFP_overlaps <- findOverlaps(GFP_peaks, promoters)
  GFP_genes <- unique(names(genes)[subjectHits(GFP_overlaps)])
  cat("Found", length(GFP_genes), "genes with GFP peaks in promoter regions\n")
  
  cat("Finding genes with SOX peaks in promoter regions\n")
  SOX_overlaps <- findOverlaps(SOX_peaks, promoters)
  SOX_genes <- unique(names(genes)[subjectHits(SOX_overlaps)])
  cat("Found", length(SOX_genes), "genes with SOX peaks in promoter regions\n")
  
  # Find genes with both GFP and SOX peaks in promoter regions
  common_genes <- intersect(GFP_genes, SOX_genes)
  cat("Found", length(common_genes), "genes with both GFP and SOX peaks in promoter regions\n")
  
  # Get gene information for common genes
  common_gene_info <- data.frame(
    gene_id = common_genes,
    gene_name = gene_names[common_genes]
  )
  
  # Sort by gene name
  common_gene_info <- common_gene_info %>% arrange(gene_name)
  
  # Save results to CSV
  write_csv(common_gene_info, output_file)
  cat("Results saved to", output_file, "\n")
}

# Run the main function (choose one implementation)
# For smaller GTF files, use the standard implementation
main()

# For very large GTF files, uncomment to use the TxDb implementation
# main_txdb()
