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
# Get the script's directory using a method compatible with Rscript
if (!is_interactive) {
  # Get the command line arguments
  args <- commandArgs(trailingOnly = FALSE)
  # Find the argument that specifies the script file path
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    # Extract the path from the --file= argument
    script_path <- sub("--file=", "", file_arg)
    # Get the directory containing the script
    script_dir <- dirname(script_path)
    # Set the working directory to the script's directory
    setwd(script_dir)
    cat("Set working directory to:", script_dir, "\n")
  } else {
    warning("Could not determine script directory. Using current working directory.")
    script_dir <- getwd()
  }
} else {
  # Fallback for interactive mode (e.g., RStudio)
  script_dir <- tryCatch(
    dirname(rstudioapi::getSourceEditorContext()$path),
    error = function(e) {
      warning("RStudio API not available. Using current working directory.")
      getwd()
    }
  )
  setwd(script_dir)
  cat("Running interactively. Set working directory to:", script_dir, "\n")
}

# Parse as I will be also running this code some time in the interactive mode 
option_list <- list(
  make_option(c("-d", "--dds_dir"), type="character", default="results/deseq2",
              help="Directory containing DESeq2 results [default: %default]"),
  make_option(c("-c", "--counts"), type="character", default="results/counts/counts.csv",
              help="Path to the original full count matrix CSV file [default: %default]"),
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
  # NOTE: The 'results_path' and 'output_path' below are example hardcoded paths.
  # Adjust 'dds_dir', 'metadata', 'output_dir', and 'counts' as needed for your interactive session.
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
    top_n_degs = 20,
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
  
  # Add default for counts file in interactive mode
  opt$counts <- "../results/counts/counts.csv" # Adjust path as needed for interactive use

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
  
  # Check if files exist
  if (!file.exists(opt$metadata)) {
    warning(paste("Metadata file does not exist:", opt$metadata))
  }
  if (!file.exists(opt$counts)) {
    warning(paste("Counts file does not exist:", opt$counts))
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

# --- Function Definitions Moved Here ---

# Function to convert gene IDs to gene symbols
get_gene_symbols <- function(gene_ids, species = opt$species) {
  # Check if we have any gene IDs to convert
  if (length(gene_ids) == 0) {
    return(character(0))
  }
  
  # Print a sample of gene IDs for debugging
  log_message(paste("Sample gene IDs:", paste(head(gene_ids, 3), collapse=", ")))
  
  # Check if these are Ensembl IDs with version numbers (e.g., ENSG00000123456.1)
  is_first_id_valid <- !is.na(gene_ids[1]) && gene_ids[1] != ""
  has_version <- FALSE
  if (is_first_id_valid) {
      has_version <- grepl("^ENS[A-Z]*[0-9]+\\.[0-9]+$", gene_ids[1])
  }

  if (has_version) {
    log_message("Detected Ensembl IDs with version numbers. Removing version numbers for mapping.")
    original_ids <- gene_ids
    gene_ids_no_version <- sub("\\.[0-9]+$", "", gene_ids)
    id_mapping <- setNames(original_ids, gene_ids_no_version)
    gene_ids_for_mapping <- gene_ids_no_version
  } else {
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
  
  # Determine primary keytype based on ID format
  primary_keytype <- "SYMBOL" # Default guess
  if (is_first_id_valid) {
      if (grepl("^ENS[A-Z]*[0-9]+", gene_ids_for_mapping[1])) {
          primary_keytype <- "ENSEMBL"
      } else if (grepl("^[NX]M_|^[NX]R_", gene_ids_for_mapping[1])) {
          primary_keytype <- "REFSEQ"
      } else if (grepl("^[0-9]+$", gene_ids_for_mapping[1])) {
           primary_keytype <- "ENTREZID"
      }
  }
  log_message(paste("Attempting mapping with primary keytype:", primary_keytype))

  # Try mapping using the determined primary keytype
  symbols <- tryCatch({
    mapIds(annotation_db, 
           keys = gene_ids_for_mapping, 
           column = "SYMBOL",
           keytype = primary_keytype, 
           multiVals = "first")
  }, error = function(e) {
    log_message(paste("Error with primary keytype", primary_keytype, ":", e$message))
    NULL
  })

  # Check success rate
  success_rate <- if (!is.null(symbols)) sum(!is.na(symbols)) / length(gene_ids_for_mapping) else 0
  log_message(paste("Mapping success rate with", primary_keytype, ":", round(success_rate * 100, 1), "%"))

  # If primary keytype failed or had low success, try alternatives
  alternative_keytypes <- c("ENSEMBL", "SYMBOL", "REFSEQ", "ENTREZID", "GENENAME")
  alternative_keytypes <- setdiff(alternative_keytypes, primary_keytype) # Remove the one we already tried

  if (success_rate < 0.5) {
      log_message("Primary keytype mapping insufficient. Trying alternatives.")
      for (key_type in alternative_keytypes) {
          log_message(paste("Trying alternative keytype:", key_type))
          symbols_alt <- tryCatch({
              mapIds(annotation_db, 
                     keys = gene_ids_for_mapping, 
                     column = "SYMBOL",
                     keytype = key_type, 
                     multiVals = "first")
          }, error = function(e) {
              log_message(paste("Error with alternative keytype", key_type, ":", e$message))
              NULL
          })

          success_rate_alt <- if (!is.null(symbols_alt)) sum(!is.na(symbols_alt)) / length(gene_ids_for_mapping) else 0
          log_message(paste("Mapping success rate with", key_type, ":", round(success_rate_alt * 100, 1), "%"))

          if (success_rate_alt > success_rate) {
              log_message(paste("Using results from alternative keytype:", key_type))
              symbols <- symbols_alt
              success_rate <- success_rate_alt
              if (success_rate > 0.8) break # Stop if we get good mapping
          }
          if (success_rate > 0.8) break # Stop if primary was good enough after all
      }
  }

  # Handle final result
  if (is.null(symbols) || success_rate == 0) {
      log_message("Could not map gene IDs to symbols. Using original IDs.")
      # Need to handle the versioned IDs correctly if they were stripped
      if (has_version) {
          return(original_ids)
      } else {
          return(gene_ids)
      }
  }

  # Map back to original IDs if versions were stripped, and fill NAs
  if (has_version) {
      final_symbols <- rep(NA_character_, length(original_ids))
      names(final_symbols) <- original_ids
      mapped_indices <- match(names(symbols), gene_ids_no_version)
      original_names_for_mapped <- original_ids[mapped_indices]
      final_symbols[original_names_for_mapped] <- symbols
      final_symbols[is.na(final_symbols)] <- names(final_symbols)[is.na(final_symbols)] # Fill NAs with original ID
      # Ensure order matches input
      final_symbols <- final_symbols[original_ids]
  } else {
      # Directly use the results, filling NAs
      final_symbols <- symbols
      na_indices <- is.na(final_symbols)
      final_symbols[na_indices] <- names(final_symbols)[na_indices]
      # Ensure order matches input
      final_symbols <- final_symbols[gene_ids]
  }

  return(final_symbols)
}


# Function to get ENTREZ IDs for pathway analysis
get_entrez_ids <- function(gene_ids, species = opt$species) {
    # Handle empty input
    if (length(gene_ids) == 0 || all(is.na(gene_ids))) {
        return(character(0))
    }

    # Select the appropriate annotation database based on species
    if (tolower(species) == "human") {
        annotation_db <- org.Hs.eg.db
    } else if (tolower(species) == "mouse") {
        annotation_db <- org.Mm.eg.db
    } else {
        warning("Unsupported species for ENTREZ ID conversion. Choose 'human' or 'mouse'.")
        return(NULL)
    }

    # Determine the most likely keytype
    first_valid_id <- na.omit(gene_ids)[1]
    key_type <- "SYMBOL" # Default guess
     if (!is.na(first_valid_id)) {
        if (grepl("^ENS[A-Z]*[0-9]+", first_valid_id)) {
            # Check for version numbers and strip if present
            if (grepl("\\.[0-9]+$", first_valid_id)) {
                 log_message("Stripping version numbers from Ensembl IDs for ENTREZ mapping.")
                 gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
            }
            key_type <- "ENSEMBL"
        } else if (grepl("^[NX]M_|^[NX]R_", first_valid_id)) {
            key_type <- "REFSEQ"
        } else if (grepl("^[0-9]+$", first_valid_id)) {
             # Already Entrez IDs? Return them directly, but ensure character type
             log_message("Input IDs look like ENTREZ IDs. Returning as is.")
             return(as.character(gene_ids))
        }
    }
    log_message(paste("Attempting ENTREZ ID mapping using keytype:", key_type))

    # Try to convert IDs to ENTREZ IDs
    entrez_ids <- tryCatch({
        mapIds(annotation_db, keys = gene_ids, column = "ENTREZID",
               keytype = key_type, multiVals = "first")
    }, error = function(e) {
        warning("Error converting gene IDs to ENTREZ IDs using keytype ", key_type, ": ", e$message)
        NULL
    })

    # If primary keytype failed, try others if appropriate
    alternative_keytypes <- c("ENSEMBL", "SYMBOL", "REFSEQ")
    alternative_keytypes <- setdiff(alternative_keytypes, key_type)

    if (is.null(entrez_ids) || sum(!is.na(entrez_ids)) < length(gene_ids) * 0.1) {
        log_message("Primary keytype mapping failed or insufficient. Trying alternatives for ENTREZ IDs.")
        for (alt_key_type in alternative_keytypes) {
            # Handle version stripping for ENSEMBL if it's the alternative
             gene_ids_alt <- gene_ids
             if (alt_key_type == "ENSEMBL" && grepl("^ENS[A-Z]*[0-9]+\\.[0-9]+$", first_valid_id)) {
                  log_message("Stripping version numbers from Ensembl IDs for alternative ENTREZ mapping.")
                  gene_ids_alt <- sub("\\.[0-9]+$", "", gene_ids)
             }

             log_message(paste("Trying alternative keytype for ENTREZ:", alt_key_type))
             entrez_ids_alt <- tryCatch({
                 mapIds(annotation_db, keys = gene_ids_alt, column = "ENTREZID",
                        keytype = alt_key_type, multiVals = "first")
             }, error = function(e) {
                 warning("Error converting gene IDs to ENTREZ IDs using alternative keytype ", alt_key_type, ": ", e$message)
                 NULL
             })

             # Use alternative if it's better
             if (!is.null(entrez_ids_alt) && sum(!is.na(entrez_ids_alt)) > sum(!is.na(entrez_ids))) {
                 log_message(paste("Using results from alternative keytype:", alt_key_type))
                 entrez_ids <- entrez_ids_alt
             }
        }
    }

    # Filter out NA values before returning
    final_entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    num_mapped <- length(final_entrez_ids)
    num_input <- length(na.omit(gene_ids))
    log_message(paste("Successfully mapped", num_mapped, "out of", num_input, "input IDs to ENTREZ IDs."))

    if (num_mapped == 0) return(NULL) # Return NULL if no IDs mapped

    return(final_entrez_ids)
}

# Function to run pathway analysis
run_pathway_analysis <- function(gene_list, species = opt$species) {
  # Ensure gene_list contains ENTREZ IDs as character strings
  gene_list <- as.character(na.omit(gene_list))

  if (is.null(gene_list) || length(gene_list) < 10) {
     log_message(paste("Skipping pathway analysis: Need at least 10 valid gene IDs, found", length(gene_list)))
    return(NULL)
  }

  # Select OrgDb based on species
  OrgDb <- NULL
  if (tolower(species) == "human") {
    OrgDb <- org.Hs.eg.db
  } else if (tolower(species) == "mouse") {
    OrgDb <- org.Mm.eg.db
  } else {
    warning("Unsupported species for pathway analysis. Choose 'human' or 'mouse'.")
    return(NULL)
  }
  
  log_message(paste("Running GO enrichment for", length(gene_list), "genes using", species, "database."))
  
  # Try to run GO enrichment
  ego <- tryCatch({
    enrichGO(gene          = gene_list,
             OrgDb         = OrgDb,
             keyType       = "ENTREZID", # Explicitly state keyType
             ont           = "BP",       # Biological Process
             pAdjustMethod = "BH",       # Benjamini-Hochberg correction
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2)        # q-value cutoff
  }, error = function(e) {
    warning("Error running GO enrichment: ", e$message)
    NULL
  })
  
  if (is.null(ego)) {
      log_message("GO enrichment resulted in NULL.")
      return(NULL)
  }

  ego_results <- as.data.frame(ego)

  if (nrow(ego_results) == 0) {
    log_message("No significant GO terms found.")
    return(NULL)
  }

  log_message(paste("Found", nrow(ego_results), "significant GO terms."))
  return(ego_results)
}

# --- End of Function Definitions ---

log_message("Starting visualization generation")

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", opt$output_dir))
}

# Define output file paths
pca_plot_file <- file.path(opt$output_dir, "pca_plot.pdf")
heatmap_file <- file.path(opt$output_dir, "sample_correlation_heatmap.pdf")
top_degs_heatmap_file <- file.path(opt$output_dir, "top_unique_degs_heatmap.pdf")
volcano_dir <- file.path(opt$output_dir, "volcano_plots")
summary_dir <- file.path(opt$output_dir, "summary_files")
comp_heatmap_dir <- file.path(opt$output_dir, "comparison_heatmaps")

# Create output directories if they don't exist
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
  log_message(paste("Created summary directory:", summary_dir))
}
if (!dir.exists(volcano_dir)) {
  dir.create(volcano_dir, recursive = TRUE)
  log_message(paste("Created volcano plots directory:", volcano_dir))
}
if (!dir.exists(comp_heatmap_dir)) {
  dir.create(comp_heatmap_dir, recursive = TRUE)
  log_message(paste("Created comparison heatmaps directory:", comp_heatmap_dir))
}

# Read full metadata
log_message(paste("Reading full metadata from:", opt$metadata))
metadata_full <- read.csv(opt$metadata, check.names = FALSE, row.names = "sample")

# Print metadata structure for debugging
log_message(paste("Metadata structure:", paste(names(metadata_full), collapse=", ")))
log_message(paste("Conditions found in metadata:", paste(unique(metadata_full$condition), collapse=", ")))
log_message(paste("Sample count by condition:"))
condition_counts <- table(metadata_full$condition)
for (cond in names(condition_counts)) {
  log_message(paste(" -", cond, ":", condition_counts[cond], "samples"))
}

# Add direct check for specific samples to avoid matching issues
log_message("Finding samples for each comparison directly from DESeqDatasets:")
for (dds_file in list.files(opt$dds_dir, pattern = "_vs_", include.dirs = TRUE, full.names = TRUE)) {
  comparison_dir <- basename(dds_file)
  if (dir.exists(dds_file) && file.exists(file.path(dds_file, "dds.rds"))) {
    dds_obj <- readRDS(file.path(dds_file, "dds.rds"))
    comparison_samples <- colnames(dds_obj)
    log_message(paste(" -", comparison_dir, ":", length(comparison_samples), "samples"))
    log_message(paste("   Samples:", paste(comparison_samples, collapse=", ")))
  }
}

# Read full count matrix
log_message(paste("Reading full count matrix from:", opt$counts))
count_matrix_full <- read.csv(opt$counts, row.names = 1, check.names = FALSE)

# Ensure samples in count matrix and metadata are aligned
samples_to_keep <- intersect(rownames(metadata_full), colnames(count_matrix_full))
if(length(samples_to_keep) == 0) {
  stop("No common samples found between metadata and count matrix.")
}
if(length(samples_to_keep) < ncol(count_matrix_full) || length(samples_to_keep) < nrow(metadata_full)) {
  log_message("Warning: Subsetting count matrix and metadata to common samples.")
}
count_matrix_full <- count_matrix_full[, samples_to_keep]
metadata_full <- metadata_full[samples_to_keep, , drop = FALSE]
metadata_full <- metadata_full %>% arrange(match(rownames(.), colnames(count_matrix_full)))

# Verify that sample order matches
if (!all(rownames(metadata_full) == colnames(count_matrix_full))) {
  stop("Sample order mismatch after filtering/ordering full metadata and counts.")
}

# Create a full DESeqDataSet for VST calculation across all samples
log_message("Creating full DESeqDataSet for VST")
dds_full <- DESeqDataSetFromMatrix(
  countData = count_matrix_full,
  colData = metadata_full,
  design = ~ condition
)

# Filter out low-count genes before VST (optional but recommended)
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ]

# Get variance stabilized transformed data for all samples
log_message("Performing variance stabilizing transformation on all samples")
vst_full <- vst(dds_full, blind = FALSE)

# Find all DESeq2 RDS files
log_message("Finding DESeq2 RDS files for comparisons")
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

# Create PCA plot using full vst data
log_message("Creating PCA plot using all samples")
pdf(pca_plot_file, width = opt$pca_width, height = opt$pca_height)

# The plotPCA function expects the intgroup to be present in colData(vst_full)
# Ensure that the 'condition' column is correctly placed in the colData
colData(vst_full)$condition <- metadata_full$condition

# Now generate the PCA data
pca_data <- plotPCA(vst_full, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = opt$pca_point_size) +
  geom_text_repel(size = opt$pca_text_size, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(opt$pca_title) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Apply selected theme
theme_func <- match.fun(paste0("theme_", opt$pca_theme))
pca_plot <- pca_plot + theme_func()

print(pca_plot)
dev.off()
log_message(paste("PCA plot saved to:", pca_plot_file))

# Create sample correlation heatmap using full vst data
log_message("Creating sample correlation heatmap using all samples")
sample_dists <- dist(t(assay(vst_full)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vst_full)
colnames(sample_dist_matrix) <- colnames(vst_full)

# Get condition colors based on full metadata
conditions_all <- unique(metadata_full$condition)
num_colors_needed <- length(conditions_all)
num_colors_available <- brewer.pal.info[opt$color_palette, "maxcolors"]
palette_colors <- brewer.pal(min(num_colors_needed, num_colors_available), opt$color_palette)
if (num_colors_needed > num_colors_available) {
  warning(paste("Number of conditions exceeds colors in palette", opt$color_palette, ". Colors will be recycled."))
  palette_colors <- rep(palette_colors, length.out = num_colors_needed)
}
condition_colors <- setNames(palette_colors, conditions_all)

# Create annotation data frame using full metadata
annotation_col <- data.frame(Condition = metadata_full$condition)
rownames(annotation_col) <- rownames(metadata_full)
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
comparison_top_degs <- list()

for (dds_file in dds_files) {
  comparison_name <- basename(dirname(dds_file))
  log_message(paste("Processing DESeq2 object:", dds_file))
  dds_obj <- readRDS(dds_file)
  
  # Get results
  res <- results(dds_obj)
  
  # Get top N DEGs by adjusted p-value for this comparison
  res_filtered <- as.data.frame(res) %>% filter(!is.na(padj))
  top_degs_comp <- rownames(res_filtered[order(res_filtered$padj), ])[1:min(opt$top_n_degs, nrow(res_filtered))]

  all_top_degs <- c(all_top_degs, list(top_degs_comp))
  comparison_top_degs[[comparison_name]] <- top_degs_comp
}

# Get unique top DEGs across all comparisons
unique_top_degs <- unique(unlist(all_top_degs))
log_message(paste("Number of unique top DEGs across all comparisons:", length(unique_top_degs)))

# Create heatmap of unique top DEGs across all comparisons
log_message("Creating heatmap of unique top DEGs using only relevant samples")

for (comparison_name in names(comparison_top_degs)) {
  top_degs_current_comp <- comparison_top_degs[[comparison_name]]
  if (length(top_degs_current_comp) > 0) {
    top_degs_current_comp_present <- intersect(top_degs_current_comp, rownames(assay(vst_full)))
    log_message(paste("Processing heatmap for:", comparison_name))
    log_message(paste("Found", length(top_degs_current_comp_present), "top DEGs for", comparison_name, "in VST data."))

    if (length(top_degs_current_comp_present) > 1) { # Need at least 2 genes to plot
      
      # Load the original DESeq2 object to get the exact samples used in this comparison
      dds_file_path <- file.path(opt$dds_dir, comparison_name, "dds.rds")
      if (file.exists(dds_file_path)) {
        dds_comparison <- readRDS(dds_file_path)
        comparison_samples <- colnames(dds_comparison)
        log_message(paste("Found", length(comparison_samples), "samples in DESeq2 object for", comparison_name))
        
        # Get the conditions from the comparison name
        conditions_in_comparison <- unlist(strsplit(comparison_name, "_vs_"))
        log_message(paste("Conditions in comparison:", paste(conditions_in_comparison, collapse=", ")))
        
        # Find samples that are both in the VST data and were used in the comparison
        relevant_samples <- intersect(comparison_samples, colnames(assay(vst_full)))
        log_message(paste("After intersection with VST data, found", length(relevant_samples), "samples"))
        
        if (length(relevant_samples) < 2) {
          # Fall back to using all samples from the conditions in the comparison
          log_message("Not enough samples found from DESeq object. Using all samples from these conditions.")
          relevant_samples <- rownames(metadata_full[metadata_full$condition %in% conditions_in_comparison, ])
          relevant_samples <- intersect(relevant_samples, colnames(assay(vst_full)))
          log_message(paste("Found", length(relevant_samples), "samples from metadata matching conditions"))
        }
        
        if (length(relevant_samples) < 2) {
          log_message(paste("Skipping heatmap for", comparison_name, ": Fewer than 2 relevant samples found."))
          next # Skip to the next comparison
        }
        
        # Print sample-condition mapping for debugging
        sample_conditions <- metadata_full[relevant_samples, "condition", drop=FALSE]
        for (i in 1:nrow(sample_conditions)) {
          log_message(paste("Sample", rownames(sample_conditions)[i], "has condition", sample_conditions[i,1]))
        }
        
        log_message(paste("Plotting heatmap for", comparison_name, "using samples:", paste(relevant_samples, collapse=", ")))
        
        # Subset the VST data for relevant genes AND relevant samples
        top_degs_counts_subset <- assay(vst_full)[top_degs_current_comp_present, relevant_samples, drop=FALSE]
        
        # Scale the subsetted counts
        top_degs_counts_scaled <- t(scale(t(top_degs_counts_subset)))
        # Handle cases where scaling might produce NaNs (e.g., zero variance)
        top_degs_counts_scaled[is.nan(top_degs_counts_scaled)] <- 0 

        # Subset the annotation data frame
        annotation_col_subset <- annotation_col[relevant_samples, , drop = FALSE]

        # Create a new annotation colors list with only the relevant conditions
        relevant_condition_colors <- ann_colors$Condition[names(ann_colors$Condition) %in% conditions_in_comparison]
        ann_colors_subset <- list(Condition = relevant_condition_colors)
        
        # Get gene symbols for rownames if possible
        heatmap_rownames <- get_gene_symbols(rownames(top_degs_counts_scaled), species = opt$species)

        # Define file path for the heatmap
        # Ensure the subdirectory exists (it should have been created in the summary loop)
        comparison_heatmap_dir_specific <- file.path(comp_heatmap_dir, comparison_name)
        if (!dir.exists(comparison_heatmap_dir_specific)) { dir.create(comparison_heatmap_dir_specific, recursive = TRUE) }
        comp_heatmap_file_pdf <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_top_unique_degs_heatmap.pdf"))
        comp_heatmap_file_png <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_top_unique_degs_heatmap.png"))

        # Generate heatmap using subsetted data and save as PDF
        pdf(comp_heatmap_file_pdf, width = opt$top_degs_width, height = opt$top_degs_height)
        pheatmap(top_degs_counts_scaled,
                 annotation_col = annotation_col_subset, # Use subsetted annotation
                 annotation_colors = ann_colors_subset,  # Use subset of annotation colors
                 labels_row = heatmap_rownames,
                 show_rownames = TRUE,
                 clustering_method = "ward.D2",
                 main = paste("Top Unique DEGs -", comparison_name), # Updated title
                 fontsize = 10,
                 fontsize_row = max(4, 10 - length(top_degs_current_comp_present) / 10))
        dev.off()
        log_message(paste("Heatmap of unique top DEGs for", comparison_name, "saved to:", comp_heatmap_file_pdf))

        # Generate heatmap using subsetted data and save as PNG
        png(comp_heatmap_file_png, width = opt$top_degs_width, height = opt$top_degs_height, units = "in", res = 300)
        pheatmap(top_degs_counts_scaled,
                 annotation_col = annotation_col_subset, # Use subsetted annotation
                 annotation_colors = ann_colors_subset,  # Use subset of annotation colors
                 labels_row = heatmap_rownames,
                 show_rownames = TRUE,
                 clustering_method = "ward.D2",
                 main = paste("Top Unique DEGs -", comparison_name), # Updated title
                 fontsize = 10,
                 fontsize_row = max(4, 10 - length(top_degs_current_comp_present) / 10))
        dev.off()
        log_message(paste("Heatmap of unique top DEGs for", comparison_name, "saved to:", comp_heatmap_file_png))

      } else {
        log_message(paste("Skipping unique top DEGs heatmap for", comparison_name, ": DESeq2 object not found."))
      }
    } else {
      log_message(paste("Skipping unique top DEGs heatmap for", comparison_name, ": Not enough unique top DEGs found in VST data."))
    }
  } else {
    log_message(paste("Skipping unique top DEGs heatmap for", comparison_name, ": No top DEGs found for this comparison."))
  }
}

# Create summary files and comparison-specific heatmaps for each comparison
log_message("Creating summary files and comparison-specific heatmaps")

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
  log_message(paste("Processing summary and heatmap for:", comparison_name))

  # Create subdirectories for this comparison
  comparison_summary_dir <- file.path(summary_dir, comparison_name)
  comparison_heatmap_dir <- file.path(comp_heatmap_dir, comparison_name)
  if (!dir.exists(comparison_summary_dir)) { dir.create(comparison_summary_dir, recursive = TRUE) }
  if (!dir.exists(comparison_heatmap_dir)) { dir.create(comparison_heatmap_dir, recursive = TRUE) }

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
  detailed_file <- file.path(comparison_summary_dir, "detailed_results.csv")
  write.csv(res_df, detailed_file, row.names = FALSE)
  log_message(paste("Detailed results saved to:", detailed_file))
  
  # Create significant genes file
  sig_genes_file <- file.path(comparison_summary_dir, "significant_genes.csv")
  sig_genes_df <- res_df[res_df$padj < opt$volcano_pval_cutoff & !is.na(res_df$padj), ]
  sig_genes_df <- sig_genes_df[order(sig_genes_df$padj), ]
  write.csv(sig_genes_df, sig_genes_file, row.names = FALSE)
  log_message(paste("Significant genes saved to:", sig_genes_file))
  
  # Create up and down regulated gene lists
  up_genes_file <- file.path(comparison_summary_dir, "up_regulated.csv")
  up_genes_df <- res_df[res_df$padj < opt$volcano_pval_cutoff & 
                        res_df$log2FoldChange > opt$volcano_fc_cutoff & 
                        !is.na(res_df$padj), ]
  up_genes_df <- up_genes_df[order(up_genes_df$log2FoldChange, decreasing = TRUE), ]
  write.csv(up_genes_df, up_genes_file, row.names = FALSE)
  log_message(paste("Up-regulated genes saved to:", up_genes_file))
  
  down_genes_file <- file.path(comparison_summary_dir, "down_regulated.csv")
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
      pathway_file <- file.path(comparison_summary_dir, "pathway_analysis.csv")
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
        up_pathway_file <- file.path(comparison_summary_dir, "up_pathway_analysis.csv")
        write.csv(up_pathway_results, up_pathway_file, row.names = FALSE)
        log_message(paste("Up-regulated pathway analysis results saved to:", up_pathway_file))
      }
    }
    
    down_entrez <- get_entrez_ids(down_genes_df$gene_id, species = opt$species)
    down_entrez <- down_entrez[!is.na(down_entrez)]
    
    if (length(down_entrez) >= 10) {
      down_pathway_results <- run_pathway_analysis(down_entrez, species = opt$species)
      
      if (!is.null(down_pathway_results) && nrow(down_pathway_results) > 0) {
        down_pathway_file <- file.path(comparison_summary_dir, "down_pathway_analysis.csv")
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
log_message("Summary file and comparison heatmap generation complete")

# Create volcano plots for each comparison
log_message("Creating volcano plots with gene symbols")

for (dds_file in dds_files) {
  comparison_name <- basename(dirname(dds_file))
  log_message(paste("Creating volcano plot for:", comparison_name))

  # Create a subdirectory for this comparison within volcano_dir
  comparison_volcano_dir <- file.path(volcano_dir, comparison_name)
  if (!dir.exists(comparison_volcano_dir)) {
    dir.create(comparison_volcano_dir, recursive = TRUE)
    log_message(paste("Created comparison volcano directory:", comparison_volcano_dir))
  }
  
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
  res_df$significant <- "Not Significant"
  res_df$significant[res_df$padj < opt$volcano_pval_cutoff & !is.na(res_df$padj) & res_df$log2FoldChange > opt$volcano_fc_cutoff] <- "Up-regulated"
  res_df$significant[res_df$padj < opt$volcano_pval_cutoff & !is.na(res_df$padj) & res_df$log2FoldChange < -opt$volcano_fc_cutoff] <- "Down-regulated"
  res_df$significant <- factor(res_df$significant, levels = c("Up-regulated", "Down-regulated", "Not Significant"))

  # Define colors for volcano plot
  volcano_colors <- c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")
  
  # Get top genes for labeling (based on p-value and fold change)
  top_genes <- res_df %>%
    filter(!is.na(padj) & padj < opt$volcano_pval_cutoff & abs(log2FoldChange) > opt$volcano_fc_cutoff) %>%
    mutate(rank_metric = -log10(padj) * abs(log2FoldChange)) %>%
    arrange(desc(rank_metric)) %>%
    head(opt$volcano_top_n)
  
  # Create volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(opt$volcano_pval_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-opt$volcano_fc_cutoff, opt$volcano_fc_cutoff), linetype = "dashed") +
    scale_color_manual(values = volcano_colors, name = "Significance") +
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
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  # Save volcano plot - with error handling
  volcano_file <- file.path(comparison_volcano_dir, "volcano_plot.pdf")
  tryCatch({
    pdf(volcano_file, width = opt$volcano_width, height = opt$volcano_height)
    print(volcano_plot)
    dev.off()
    log_message(paste("Volcano plot saved to:", volcano_file))
  }, error = function(e) {
    warning(paste("Failed to save PDF volcano plot:", e$message))
  })
  
  # Also save as PNG for easier viewing - with error handling
  volcano_png <- file.path(comparison_volcano_dir, "volcano_plot.png")
  tryCatch({
    ggsave(volcano_png, volcano_plot, width = opt$volcano_width, height = opt$volcano_height, dpi = 300)
    log_message(paste("Volcano plot saved to:", volcano_png))
  }, error = function(e) {
    warning(paste("Failed to save PNG volcano plot:", e$message))
  })
}

log_message("Visualization generation complete")
