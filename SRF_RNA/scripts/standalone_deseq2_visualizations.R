#!/usr/bin/env Rscript

# This script generates visualizations from DESeq2 results
# It creates PCA plots, sample correlation heatmaps, and heatmaps of top DEGs
# This is a standalone version that can be run independently of the Snakefile

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
  library(org.Hs.eg.db)  # For human gene annotation
  library(org.Mm.eg.db)  # For mouse gene annotation
  library(clusterProfiler) # For pathway analysis
  library(DOSE)          # For disease ontology
})

# Check if script is being run interactively
is_interactive <- interactive()

# Set working directory to script's location
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_dir)

# Parse as I will be also running this code some time in the interactive mode 
option_list <- list(
  make_option(c("-d", "--dds_dir"), type="character", default="results/deseq2",
              help="Directory containing DESeq2 results [default: %default]"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",
              help="Metadata file [default: %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="results/deseq2",
              help="Output directory for visualization files [default: %default]"),
  make_option(c("--pca_width"), type="numeric", default=10,
              help="Width of PCA plot in inches [default: %default]"),
  make_option(c("--pca_height"), type="numeric", default=8,
              help="Height of PCA plot in inches [default: %default]"),
  make_option(c("--heatmap_width"), type="numeric", default=10,
              help="Width of sample correlation heatmap in inches [default: %default]"),
  make_option(c("--heatmap_height"), type="numeric", default=8,
              help="Height of sample correlation heatmap in inches [default: %default]"),
  make_option(c("--top_degs_width"), type="numeric", default=12,
              help="Width of top DEGs heatmap in inches [default: %default]"),
  make_option(c("--top_degs_height"), type="numeric", default=14,
              help="Height of top DEGs heatmap in inches [default: %default]"),
  make_option(c("--top_n_degs"), type="numeric", default=50,
              help="Number of top DEGs to include in heatmap [default: %default]"),
  make_option(c("--pca_title"), type="character", default="PCA of RNA-seq Samples",
              help="Title for PCA plot [default: %default]"),
  make_option(c("--heatmap_title"), type="character", default="Sample Correlation Heatmap",
              help="Title for sample correlation heatmap [default: %default]"),
  make_option(c("--top_degs_title"), type="character", default="Top Differentially Expressed Genes",
              help="Title for top DEGs heatmap [default: %default]"),
  make_option(c("--pca_point_size"), type="numeric", default=3,
              help="Point size for PCA plot [default: %default]"),
  make_option(c("--pca_text_size"), type="numeric", default=4,
              help="Text size for PCA plot labels [default: %default]"),
  make_option(c("--pca_theme"), type="character", default="bw",
              help="Theme for PCA plot (bw, classic, minimal, etc.) [default: %default]"),
  make_option(c("--color_palette"), type="character", default="Set1",
              help="RColorBrewer palette for condition colors [default: %default]"),
  make_option(c("--volcano_width"), type="numeric", default=10,
              help="Width of volcano plot in inches [default: %default]"),
  make_option(c("--volcano_height"), type="numeric", default=8,
              help="Height of volcano plot in inches [default: %default]"),
  make_option(c("--volcano_top_n"), type="numeric", default=20,
              help="Number of top genes to label in volcano plot [default: %default]"),
  make_option(c("--volcano_fc_cutoff"), type="numeric", default=1,
              help="Log2 fold change cutoff for volcano plot [default: %default]"),
  make_option(c("--volcano_pval_cutoff"), type="numeric", default=0.05,
              help="Adjusted p-value cutoff for volcano plot [default: %default]"),
  make_option(c("--species"), type="character", default="human",
              help="Species for gene annotation (human or mouse) [default: %default]"),
  make_option(c("--run_pathway"), type="logical", default=FALSE,
              help="Run pathway analysis [default: %default]"),
  make_option(c("--verbose"), type="logical", default=TRUE,
              help="Print verbose output [default: %default]")
)

results_path <- "E:/SRF_DATA/Antonella/RNA/results/deseq2"
output_path <- "E:/SRF_CODE/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2"

# Handle command line arguments differently based on interactive mode
if (is_interactive) {
  # In interactive mode, set up default values that can be manually modified
  opt <- list(
    dds_dir = results_path,
    metadata = "../metadata.csv",
    output_dir = output_path,
    pca_width = 10,
    pca_height = 8,
    heatmap_width = 10,
    heatmap_height = 8,
    top_degs_width = 12,
    top_degs_height = 14,
    top_n_degs = 50,
    pca_title = "PCA of RNA-seq Samples",
    heatmap_title = "Sample Correlation Heatmap",
    top_degs_title = "Top Differentially Expressed Genes",
    pca_point_size = 3,
    pca_text_size = 4,
    pca_theme = "bw",
    color_palette = "Set1",
    volcano_width = 10,
    volcano_height = 8,
    volcano_top_n = 20,
    volcano_fc_cutoff = 1,
    volcano_pval_cutoff = 0.05,
    species = "human",
    run_pathway = FALSE,
    verbose = TRUE
  )
  
  # Print message about interactive mode
  cat("Running in interactive mode. You can modify the 'opt' list before proceeding.\n")
  cat("Example: opt$dds_dir <- 'path/to/deseq2/results'\n")
  cat("         opt$metadata <- 'path/to/metadata.csv'\n")
  cat("         opt$top_n_degs <- 100\n")
  cat("         opt$species <- 'mouse'  # For mouse gene annotation\n")
  cat("         opt$run_pathway <- TRUE # Run pathway analysis\n")
  cat("Continue with your analysis after setting your desired parameters.\n\n")
} else {
  # In non-interactive mode, parse command line arguments
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
}

# Function to print messages if verbose is TRUE
log_message <- function(message) {
  if (opt$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  }
}

# Add a function to validate parameters
validate_parameters <- function() {
  # Check if directories exist
  if (!dir.exists(opt$dds_dir)) {
    warning(paste("DESeq2 directory does not exist:", opt$dds_dir))
  }
  
  # Check if metadata file exists
  if (!file.exists(opt$metadata)) {
    warning(paste("Metadata file does not exist:", opt$metadata))
  }
  
  # Check numeric parameters are positive
  numeric_params <- c("pca_width", "pca_height", "heatmap_width", "heatmap_height", 
                     "top_degs_width", "top_degs_height", "top_n_degs", 
                     "pca_point_size", "pca_text_size")
  
  for (param in numeric_params) {
    if (opt[[param]] <= 0) {
      warning(paste(param, "should be positive, current value:", opt[[param]]))
    }
  }
  
  # Check if color palette exists
  if (!opt$color_palette %in% rownames(brewer.pal.info)) {
    warning(paste("Color palette not found:", opt$color_palette, 
                 "Available palettes:", paste(rownames(brewer.pal.info), collapse=", ")))
  }
  
  # Check if PCA theme is valid
  valid_themes <- c("bw", "classic", "minimal", "light", "dark", "gray")
  if (!opt$pca_theme %in% valid_themes) {
    warning(paste("PCA theme may not be valid:", opt$pca_theme,
                 "Recommended themes:", paste(valid_themes, collapse=", ")))
  }
}

# Validate parameters
if (is_interactive) {
  cat("You can run validate_parameters() to check your parameter settings before proceeding.\n")
} else {
  validate_parameters()
}

log_message("Starting visualization generation")

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", opt$output_dir))
}

# Define output file paths
pca_plot_file <- file.path(opt$output_dir, "pca_plot.pdf")
heatmap_file <- file.path(opt$output_dir, "sample_correlation_heatmap.pdf")
top_degs_heatmap_file <- file.path(opt$output_dir, "top_degs_heatmap.pdf")
volcano_dir <- file.path(opt$output_dir, "volcano_plots")
summary_dir <- file.path(opt$output_dir, "summary_files")

# Create summary directory if it doesn't exist
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
  log_message(paste("Created summary directory:", summary_dir))
}

# Read metadata
log_message(paste("Reading metadata from:", opt$metadata))
metadata <- read.csv(opt$metadata, check.names = FALSE)

# Find all DESeq2 RDS files
log_message("Finding DESeq2 RDS files")
comparisons <- list.dirs(opt$dds_dir, full.names = TRUE, recursive = FALSE)
comparisons <- comparisons[grepl("_vs_", basename(comparisons))]
dds_files <- file.path(comparisons, "dds.rds")
dds_files <- dds_files[file.exists(dds_files)]

if (length(dds_files) == 0) {
  stop("No DESeq2 RDS files found in ", opt$dds_dir)
}

log_message(paste("Found", length(dds_files), "DESeq2 RDS files"))
for (file in dds_files) {
  log_message(paste("  -", file))
}

# Load the first DESeq2 object to get normalized counts for all samples
log_message(paste("Loading DESeq2 object from:", dds_files[1]))
dds <- readRDS(dds_files[1])

# Get variance stabilized transformed data for all samples
log_message("Performing variance stabilizing transformation")
vst <- vst(dds, blind = FALSE)

# Create PCA plot
log_message("Creating PCA plot")
pdf(pca_plot_file, width = opt$pca_width, height = opt$pca_height)
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = opt$pca_point_size) +
  geom_text_repel(size = opt$pca_text_size, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(opt$pca_title) +
  theme_get() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Apply selected theme
if (opt$pca_theme == "bw") {
  pca_plot <- pca_plot + theme_bw()
} else if (opt$pca_theme == "classic") {
  pca_plot <- pca_plot + theme_classic()
} else if (opt$pca_theme == "minimal") {
  pca_plot <- pca_plot + theme_minimal()
}

print(pca_plot)
dev.off()
log_message(paste("PCA plot saved to:", pca_plot_file))

# Create sample correlation heatmap
log_message("Creating sample correlation heatmap")
sample_dists <- dist(t(assay(vst)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst)
colnames(sample_dist_matrix) <- colnames(vst)

# Get condition colors
condition_colors <- brewer.pal(min(length(unique(metadata$condition)), 9), opt$color_palette)
names(condition_colors) <- unique(metadata$condition)
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- metadata$sample
ann_colors <- list(Condition = condition_colors)

pdf(heatmap_file, width = opt$heatmap_width, height = opt$heatmap_height)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = opt$heatmap_title,
         fontsize = 10)
dev.off()
log_message(paste("Sample correlation heatmap saved to:", heatmap_file))

# Collect top DEGs from all comparisons
log_message("Collecting top DEGs from all comparisons")
all_top_degs <- list()

for (dds_file in dds_files) {
  log_message(paste("Processing DESeq2 object:", dds_file))
  dds_obj <- readRDS(dds_file)
  
  # Get results
  res <- results(dds_obj)
  
  # Get top N DEGs by adjusted p-value
  top_degs <- rownames(res[order(res$padj), ])[1:opt$top_n_degs]
  all_top_degs <- c(all_top_degs, list(top_degs))
}

# Get unique top DEGs
unique_top_degs <- unique(unlist(all_top_degs))
log_message(paste("Number of unique top DEGs:", length(unique_top_degs)))

# Create heatmap of top DEGs
log_message("Creating heatmap of top DEGs")
# Extract normalized counts for top DEGs
top_degs_counts <- assay(vst)[unique_top_degs, ]

# Scale the counts for better visualization
top_degs_counts_scaled <- t(scale(t(top_degs_counts)))

pdf(top_degs_heatmap_file, width = opt$top_degs_width, height = opt$top_degs_height)
pheatmap(top_degs_counts_scaled,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         clustering_method = "ward.D2",
         main = opt$top_degs_title,
         fontsize = 10)
dev.off()
log_message(paste("Top DEGs heatmap saved to:", top_degs_heatmap_file))

# Function to convert gene IDs to gene symbols
get_gene_symbols <- function(gene_ids, species = opt$species) {
  # Check if we have any gene IDs to convert
  if (length(gene_ids) == 0) {
    return(character(0))
  }
  
  # Print a sample of gene IDs for debugging
  log_message(paste("Sample gene IDs:", paste(head(gene_ids, 3), collapse=", ")))
  
  # Check if these are Ensembl IDs with version numbers (e.g., ENSG00000123456.1)
  has_version <- grepl("^ENS[A-Z]*[0-9]+\\.[0-9]+$", gene_ids[1])
  
  if (has_version) {
    log_message("Detected Ensembl IDs with version numbers. Removing version numbers for mapping.")
    # Create a mapping to keep track of original IDs
    original_ids <- gene_ids
    # Remove version numbers (everything after the dot)
    gene_ids_no_version <- sub("\\.[0-9]+$", "", gene_ids)
    # Create a mapping from stripped IDs back to original IDs
    id_mapping <- setNames(original_ids, gene_ids_no_version)
    # Use the stripped IDs for mapping
    gene_ids_for_mapping <- gene_ids_no_version
  } else {
    # Use original IDs for mapping
    gene_ids_for_mapping <- gene_ids
  }
  
  # Select the appropriate annotation database based on species
  if (tolower(species) == "human") {
    annotation_db <- org.Hs.eg.db
  } else if (tolower(species) == "mouse") {
    annotation_db <- org.Mm.eg.db
  } else {
    warning("Unsupported species. Choose 'human' or 'mouse'.")
    return(gene_ids)  # Return original IDs
  }
  
  # Try to map using ENSEMBL keytype first (most common for RNA-seq data)
  log_message("Trying ENSEMBL keytype for gene ID mapping")
  symbols <- tryCatch({
    mapped_symbols <- mapIds(annotation_db, 
                            keys = gene_ids_for_mapping, 
                            column = "SYMBOL",
                            keytype = "ENSEMBL", 
                            multiVals = "first")
    
    # Check if mapping was successful
    if (sum(!is.na(mapped_symbols)) > length(gene_ids_for_mapping) * 0.1) {
      log_message(paste("Successfully mapped", sum(!is.na(mapped_symbols)), "out of", length(gene_ids_for_mapping), "gene IDs"))
      mapped_symbols
    } else {
      log_message("Less than 10% of IDs mapped. Trying other keytypes.")
      NULL
    }
  }, error = function(e) {
    log_message(paste("Error with ENSEMBL keytype:", e$message))
    NULL
  })
  
  # If ENSEMBL mapping failed, try other keytypes
  if (is.null(symbols)) {
    key_types <- c("SYMBOL", "REFSEQ", "ENTREZID", "GENENAME")
    
    for (key_type in key_types) {
      log_message(paste("Trying", key_type, "keytype for gene ID mapping"))
      
      symbols <- tryCatch({
        mapped_symbols <- mapIds(annotation_db, 
                                keys = gene_ids_for_mapping, 
                                column = "SYMBOL",
                                keytype = key_type, 
                                multiVals = "first")
        
        # Check if mapping was successful
        if (sum(!is.na(mapped_symbols)) > length(gene_ids_for_mapping) * 0.1) {
          log_message(paste("Successfully mapped", sum(!is.na(mapped_symbols)), "out of", length(gene_ids_for_mapping), "gene IDs"))
          mapped_symbols
          break
        } else {
          log_message("Less than 10% of IDs mapped. Trying next keytype.")
          NULL
        }
      }, error = function(e) {
        log_message(paste("Error with", key_type, "keytype:", e$message))
        NULL
      })
      
      if (!is.null(symbols)) {
        break
      }
    }
  }
  
  # If we still don't have symbols, return original IDs
  if (is.null(symbols)) {
    log_message("Could not map gene IDs to symbols. Using original IDs.")
    return(gene_ids)
  }
  
  # Replace NAs with original IDs
  if (has_version) {
    # For IDs with version numbers, we need to map back to the original IDs
    names(symbols) <- id_mapping[names(symbols)]
    result <- symbols
    result[is.na(result)] <- names(result)[is.na(result)]
  } else {
    result <- symbols
    result[is.na(result)] <- names(result)[is.na(result)]
  }
  
  return(result)
}

# Function to get ENTREZ IDs for pathway analysis
get_entrez_ids <- function(gene_ids, species = opt$species) {
  # Select the appropriate annotation database based on species
  if (tolower(species) == "human") {
    annotation_db <- org.Hs.eg.db
    key_type <- "ENSEMBL"
    if (!grepl("^ENSG", gene_ids[1])) {
      # If not ENSEMBL IDs, try to determine the ID type
      if (grepl("^NM_|^XM_|^NR_|^XR_", gene_ids[1])) {
        key_type <- "REFSEQ"
      } else {
        key_type <- "SYMBOL"  # Default to symbol if can't determine
      }
    }
  } else if (tolower(species) == "mouse") {
    annotation_db <- org.Mm.eg.db
    key_type <- "ENSEMBL"
    if (!grepl("^ENSMUS", gene_ids[1])) {
      # If not ENSEMBL IDs, try to determine the ID type
      if (grepl("^NM_|^XM_|^NR_|^XR_", gene_ids[1])) {
        key_type <- "REFSEQ"
      } else {
        key_type <- "SYMBOL"  # Default to symbol if can't determine
      }
    }
  } else {
    stop("Unsupported species. Choose 'human' or 'mouse'.")
  }
  
  # Try to convert IDs to ENTREZ IDs
  tryCatch({
    entrez_ids <- mapIds(annotation_db, keys = gene_ids, column = "ENTREZID", 
                        keytype = key_type, multiVals = "first")
    return(entrez_ids)
  }, error = function(e) {
    warning("Error converting gene IDs to ENTREZ IDs: ", e$message)
    return(NULL)  # Return NULL if conversion fails
  })
}

# Function to run pathway analysis
run_pathway_analysis <- function(gene_list, species = opt$species) {
  if (is.null(gene_list) || length(gene_list) < 10) {
    return(NULL)
  }
  
  # Try to run GO enrichment
  tryCatch({
    if (tolower(species) == "human") {
      ego <- enrichGO(gene = gene_list,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
    } else {
      ego <- enrichGO(gene = gene_list,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
    }
    
    if (is.null(ego) || nrow(ego) == 0) {
      return(NULL)
    }
    
    return(as.data.frame(ego))
  }, error = function(e) {
    warning("Error running GO enrichment: ", e$message)
    return(NULL)
  })
}

# Create summary files for each comparison
log_message("Creating summary files with DESeq2 results")

# Create a master summary file
master_summary <- data.frame(
  comparison = character(),
  total_genes = integer(),
  up_regulated = integer(),
  down_regulated = integer(),
  significant_genes = integer(),
  max_log2fc = numeric(),
  min_log2fc = numeric(),
  top_up_gene = character(),
  top_down_gene = character(),
  stringsAsFactors = FALSE
)

for (dds_file in dds_files) {
  comparison_name <- basename(dirname(dds_file))
  log_message(paste("Creating summary for:", comparison_name))
  
  # Load DESeq2 results
  dds_obj <- readRDS(dds_file)
  res <- results(dds_obj)
  
  # Convert to data frame
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # Add gene symbols
  log_message("Adding gene symbols to results")
  res_df$gene_symbol <- get_gene_symbols(res_df$gene_id, species = opt$species)
  
  # Calculate summary statistics
  total_genes <- nrow(res_df)
  significant_genes <- sum(res_df$padj < opt$volcano_pval_cutoff & !is.na(res_df$padj))
  up_regulated <- sum(res_df$padj < opt$volcano_pval_cutoff & res_df$log2FoldChange > opt$volcano_fc_cutoff & !is.na(res_df$padj))
  down_regulated <- sum(res_df$padj < opt$volcano_pval_cutoff & res_df$log2FoldChange < -opt$volcano_fc_cutoff & !is.na(res_df$padj))
  
  # Get top up and down regulated genes
  sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < opt$volcano_pval_cutoff, ]
  
  if (nrow(sig_genes) > 0) {
    top_up_gene <- sig_genes[which.max(sig_genes$log2FoldChange), "gene_symbol"]
    top_down_gene <- sig_genes[which.min(sig_genes$log2FoldChange), "gene_symbol"]
    max_log2fc <- max(sig_genes$log2FoldChange)
    min_log2fc <- min(sig_genes$log2FoldChange)
  } else {
    top_up_gene <- "None"
    top_down_gene <- "None"
    max_log2fc <- NA
    min_log2fc <- NA
  }
  
  # Add to master summary
  master_summary <- rbind(master_summary, data.frame(
    comparison = comparison_name,
    total_genes = total_genes,
    up_regulated = up_regulated,
    down_regulated = down_regulated,
    significant_genes = significant_genes,
    max_log2fc = max_log2fc,
    min_log2fc = min_log2fc,
    top_up_gene = top_up_gene,
    top_down_gene = top_down_gene,
    stringsAsFactors = FALSE
  ))
  
  # Create detailed results file with gene symbols
  detailed_file <- file.path(summary_dir, paste0(comparison_name, "_detailed_results.csv"))
  write.csv(res_df, detailed_file, row.names = FALSE)
  log_message(paste("Detailed results saved to:", detailed_file))
  
  # Create significant genes file
  sig_genes_file <- file.path(summary_dir, paste0(comparison_name, "_significant_genes.csv"))
  sig_genes_df <- res_df[res_df$padj < opt$volcano_pval_cutoff & !is.na(res_df$padj), ]
  sig_genes_df <- sig_genes_df[order(sig_genes_df$padj), ]
  write.csv(sig_genes_df, sig_genes_file, row.names = FALSE)
  log_message(paste("Significant genes saved to:", sig_genes_file))
  
  # Create up and down regulated gene lists
  up_genes_file <- file.path(summary_dir, paste0(comparison_name, "_up_regulated.csv"))
  up_genes_df <- res_df[res_df$padj < opt$volcano_pval_cutoff & 
                        res_df$log2FoldChange > opt$volcano_fc_cutoff & 
                        !is.na(res_df$padj), ]
  up_genes_df <- up_genes_df[order(up_genes_df$log2FoldChange, decreasing = TRUE), ]
  write.csv(up_genes_df, up_genes_file, row.names = FALSE)
  log_message(paste("Up-regulated genes saved to:", up_genes_file))
  
  down_genes_file <- file.path(summary_dir, paste0(comparison_name, "_down_regulated.csv"))
  down_genes_df <- res_df[res_df$padj < opt$volcano_pval_cutoff & 
                          res_df$log2FoldChange < -opt$volcano_fc_cutoff & 
                          !is.na(res_df$padj), ]
  down_genes_df <- down_genes_df[order(down_genes_df$log2FoldChange), ]
  write.csv(down_genes_df, down_genes_file, row.names = FALSE)
  log_message(paste("Down-regulated genes saved to:", down_genes_file))
  
  # Run pathway analysis if requested
  if (opt$run_pathway) {
    log_message("Running pathway analysis")
    
    # Get ENTREZ IDs for significant genes
    sig_entrez <- get_entrez_ids(sig_genes_df$gene_id, species = opt$species)
    sig_entrez <- sig_entrez[!is.na(sig_entrez)]
    
    # Run pathway analysis
    pathway_results <- run_pathway_analysis(sig_entrez, species = opt$species)
    
    if (!is.null(pathway_results) && nrow(pathway_results) > 0) {
      pathway_file <- file.path(summary_dir, paste0(comparison_name, "_pathway_analysis.csv"))
      write.csv(pathway_results, pathway_file, row.names = FALSE)
      log_message(paste("Pathway analysis results saved to:", pathway_file))
    } else {
      log_message("No significant pathway enrichment found")
    }
    
    # Separate pathway analysis for up and down regulated genes
    up_entrez <- get_entrez_ids(up_genes_df$gene_id, species = opt$species)
    up_entrez <- up_entrez[!is.na(up_entrez)]
    
    if (length(up_entrez) >= 10) {
      up_pathway_results <- run_pathway_analysis(up_entrez, species = opt$species)
      
      if (!is.null(up_pathway_results) && nrow(up_pathway_results) > 0) {
        up_pathway_file <- file.path(summary_dir, paste0(comparison_name, "_up_pathway_analysis.csv"))
        write.csv(up_pathway_results, up_pathway_file, row.names = FALSE)
        log_message(paste("Up-regulated pathway analysis results saved to:", up_pathway_file))
      }
    }
    
    down_entrez <- get_entrez_ids(down_genes_df$gene_id, species = opt$species)
    down_entrez <- down_entrez[!is.na(down_entrez)]
    
    if (length(down_entrez) >= 10) {
      down_pathway_results <- run_pathway_analysis(down_entrez, species = opt$species)
      
      if (!is.null(down_pathway_results) && nrow(down_pathway_results) > 0) {
        down_pathway_file <- file.path(summary_dir, paste0(comparison_name, "_down_pathway_analysis.csv"))
        write.csv(down_pathway_results, down_pathway_file, row.names = FALSE)
        log_message(paste("Down-regulated pathway analysis results saved to:", down_pathway_file))
      }
    }
  }
}

# Write master summary file
master_summary_file <- file.path(summary_dir, "master_summary.csv")
write.csv(master_summary, master_summary_file, row.names = FALSE)
log_message(paste("Master summary saved to:", master_summary_file))

# Create a README file with descriptions of all output files
readme_file <- file.path(summary_dir, "README.txt")
readme_content <- c(
  "DESeq2 Analysis Summary Files",
  "===========================",
  "",
  "This directory contains summary files from the DESeq2 differential expression analysis.",
  "",
  "File Descriptions:",
  "----------------",
  "",
  "master_summary.csv - Overview of all comparisons with key statistics",
  "  - comparison: Name of the comparison",
  "  - total_genes: Total number of genes analyzed",
  "  - up_regulated: Number of significantly up-regulated genes",
  "  - down_regulated: Number of significantly down-regulated genes",
  "  - significant_genes: Total number of significantly differentially expressed genes",
  "  - max_log2fc: Maximum log2 fold change among significant genes",
  "  - min_log2fc: Minimum log2 fold change among significant genes",
  "  - top_up_gene: Gene symbol with highest positive fold change",
  "  - top_down_gene: Gene symbol with highest negative fold change",
  "",
  "For each comparison, the following files are generated:",
  "",
  "1. [comparison]_detailed_results.csv - Complete DESeq2 results for all genes",
  "   - gene_id: Original gene identifier",
  "   - gene_symbol: Corresponding gene symbol",
  "   - baseMean: Average normalized count across all samples",
  "   - log2FoldChange: Log2 fold change",
  "   - lfcSE: Standard error of log2 fold change",
  "   - stat: Wald statistic",
  "   - pvalue: Raw p-value",
  "   - padj: Adjusted p-value (FDR)",
  "",
  "2. [comparison]_significant_genes.csv - Genes with adjusted p-value < cutoff",
  "",
  "3. [comparison]_up_regulated.csv - Significantly up-regulated genes",
  "",
  "4. [comparison]_down_regulated.csv - Significantly down-regulated genes",
  "",
  if (opt$run_pathway) c(
    "5. [comparison]_pathway_analysis.csv - Pathway enrichment analysis for all significant genes",
    "",
    "6. [comparison]_up_pathway_analysis.csv - Pathway enrichment for up-regulated genes",
    "",
    "7. [comparison]_down_pathway_analysis.csv - Pathway enrichment for down-regulated genes",
    ""
  ) else character(0),
  paste("Analysis performed with cutoffs: adjusted p-value <", opt$volcano_pval_cutoff, 
        "and absolute log2 fold change >", opt$volcano_fc_cutoff),
  paste("Date:", Sys.Date())
)

writeLines(readme_content, readme_file)
log_message(paste("README file saved to:", readme_file))

log_message("Summary file generation complete")

# Create volcano plots directory if it doesn't exist
if (!dir.exists(volcano_dir)) {
  dir.create(volcano_dir, recursive = TRUE)
  log_message(paste("Created volcano plots directory:", volcano_dir))
}

# Create summary directory if it doesn't exist
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
  log_message(paste("Created summary directory:", summary_dir))
}

# Function to check if a directory is writable
is_writable <- function(dir_path) {
  test_file <- file.path(dir_path, "test_write_permission.tmp")
  result <- tryCatch({
    file.create(test_file)
  }, error = function(e) {
    return(FALSE)
  })
  
  if (file.exists(test_file)) {
    file.remove(test_file)
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Check if output directories are writable
if (!is_writable(opt$output_dir)) {
  warning(paste("Output directory is not writable:", opt$output_dir))
  warning("Will try to use a temporary directory instead")
  temp_dir <- tempdir()
  opt$output_dir <- file.path(temp_dir, "deseq2_results")
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Update output file paths
  pca_plot_file <- file.path(opt$output_dir, "pca_plot.pdf")
  heatmap_file <- file.path(opt$output_dir, "sample_correlation_heatmap.pdf")
  top_degs_heatmap_file <- file.path(opt$output_dir, "top_degs_heatmap.pdf")
  volcano_dir <- file.path(opt$output_dir, "volcano_plots")
  summary_dir <- file.path(opt$output_dir, "summary_files")
  
  # Create directories
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_message(paste("Using temporary directory instead:", opt$output_dir))
}

# Create volcano plots for each comparison
log_message("Creating volcano plots with gene symbols")

for (dds_file in dds_files) {
  comparison_name <- basename(dirname(dds_file))
  log_message(paste("Creating volcano plot for:", comparison_name))
  
  # Load DESeq2 results
  dds_obj <- readRDS(dds_file)
  res <- results(dds_obj)
  
  # Convert to data frame for ggplot
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  # Add gene symbols
  log_message("Converting gene IDs to gene symbols")
  res_df$gene_symbol <- get_gene_symbols(res_df$gene_id, species = opt$species)
  
  # Add significance column
  res_df$significant <- ifelse(res_df$padj < opt$volcano_pval_cutoff & 
                              abs(res_df$log2FoldChange) > opt$volcano_fc_cutoff, 
                              "Significant", "Not Significant")
  
  # Get top genes for labeling
  top_genes <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(opt$volcano_top_n)
  
  # Create volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    geom_hline(yintercept = -log10(opt$volcano_pval_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-opt$volcano_fc_cutoff, opt$volcano_fc_cutoff), linetype = "dashed") +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    geom_text_repel(data = top_genes, 
                   aes(label = gene_symbol),
                   box.padding = 0.5,
                   max.overlaps = 20,
                   size = 3) +
    labs(title = paste("Volcano Plot -", comparison_name),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom")
  
  # Save volcano plot - with error handling
  volcano_file <- file.path(volcano_dir, paste0(comparison_name, "_volcano.pdf"))
  tryCatch({
    pdf(volcano_file, width = opt$volcano_width, height = opt$volcano_height)
    print(volcano_plot)
    dev.off()
    log_message(paste("Volcano plot saved to:", volcano_file))
  }, error = function(e) {
    warning(paste("Failed to save PDF volcano plot:", e$message))
  })
  
  # Also save as PNG for easier viewing - with error handling
  volcano_png <- file.path(volcano_dir, paste0(comparison_name, "_volcano.png"))
  tryCatch({
    ggsave(volcano_png, volcano_plot, width = opt$volcano_width, height = opt$volcano_height, dpi = 300)
    log_message(paste("Volcano plot saved to:", volcano_png))
  }, error = function(e) {
    warning(paste("Failed to save PNG volcano plot:", e$message))
  })
}

log_message("Visualization generation complete")
