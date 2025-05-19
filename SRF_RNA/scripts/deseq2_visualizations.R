#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)  # For human gene annotation
  library(clusterProfiler) # For pathway analysis
  library(DOSE)          # For disease ontology
})

# Define default parameters
# These can be overridden by command line arguments if the script is run non-interactively
PARAMS <- list(
  # Path to the directory containing DESeq2 comparison subdirectories (e.g., YAF_vs_GFP)
  dds_dir = "SRF_RNA/results/deseq2",
  # Path to the metadata file relative to the workspace root
  metadata = "SRF_RNA/metadata.csv",
  # Base output directory for visualizations relative to the workspace root
  output_dir = "SRF_RNA/results/viz_deseq2",
  # Path to the full count matrix file relative to the workspace root
  counts = "SRF_RNA/results/counts/all_samples_counts.txt",
  pca_width = 10,
  pca_height = 8,
  heatmap_width = 10,
  heatmap_height = 8,
  top_degs_width = 12,
  top_degs_height = 14,
  top_n_degs = 50, # Changed default from 20 to 50
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
  species = "human", # Supported: "human", "mouse"
  run_pathway = FALSE, # Set to TRUE to run GO enrichment analysis
  verbose = TRUE # Set to FALSE to suppress log messages
)

# Note: Command line argument parsing would typically go here in a non-interactive script
# For this cleanup, we assume parameters are set in the list above or handled externally.

# Function to print messages if verbose is TRUE
log_message <- function(message) {
  if (PARAMS$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  }
}

# Add a function to validate parameters
validate_parameters <- function() {
  log_message("Validating parameters...")

  # Check if base directories exist
  if (!dir.exists(PARAMS$dds_dir)) {
    warning(paste("Base DESeq2 results directory does not exist:", PARAMS$dds_dir))
  }
  if (!dir.exists(PARAMS$output_dir)) {
     log_message(paste("Output base directory does not exist:", PARAMS$output_dir, "- will attempt to create."))
  }

  # Check if files exist
  if (!file.exists(PARAMS$metadata)) {
    warning(paste("Metadata file does not exist:", PARAMS$metadata))
  }
  if (!file.exists(PARAMS$counts)) {
    warning(paste("Counts file does not exist:", PARAMS$counts))
  }

  # Check numeric parameters are positive
  numeric_params <- c("pca_width", "pca_height", "heatmap_width", "heatmap_height",
                     "top_degs_width", "top_degs_height", "top_n_degs",
                     "pca_point_size", "pca_text_size", "volcano_width", "volcano_height", "volcano_top_n", "volcano_fc_cutoff", "volcano_pval_cutoff")

  for (param in numeric_params) {
    if (!is.numeric(PARAMS[[param]]) || PARAMS[[param]] <= 0) {
      warning(paste(param, "should be a positive number, current value:", PARAMS[[param]]))
    }
  }

  # Check if color palette exists
  if (!PARAMS$color_palette %in% rownames(brewer.pal.info)) {
    warning(paste("Color palette not found:", PARAMS$color_palette,
                 "Available palettes:", paste(rownames(brewer.pal.info), collapse=", ")))
  }

  # Check if PCA theme is valid
  valid_themes <- c("bw", "classic", "minimal", "light", "dark", "gray")
  if (!PARAMS$pca_theme %in% valid_themes) {
    warning(paste("PCA theme may not be valid:", PARAMS$pca_theme,
                 "Recommended themes:", paste(valid_themes, collapse=", ")))
  }

  # Check species
  if (!tolower(PARAMS$species) %in% c("human", "mouse")) {
       warning(paste("Unsupported species:", PARAMS$species, ". Supported species are 'human' and 'mouse'."))
  }

  log_message("Parameter validation complete.")
}

# Validate parameters if not in interactive mode
if (!interactive()) {
  validate_parameters()
} else {
   cat("Note: Running in interactive mode. Parameters are defined in the PARAMS list.\n")
   cat("You can run validate_parameters() to check your parameter settings before proceeding.\n")
}

# --- Function Definitions ---

# Function to convert gene IDs to gene symbols
get_gene_symbols <- function(gene_ids, species = PARAMS$species) {
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
  
  annotation_db <- org.Hs.eg.db
  
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
get_entrez_ids <- function(gene_ids, species = PARAMS$species) {
    # Handle empty input
    if (length(gene_ids) == 0 || all(is.na(gene_ids))) {
        return(character(0))
    }

    annotation_db <- org.Hs.eg.db

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
run_pathway_analysis <- function(gene_list, species = PARAMS$species) {
  # Ensure gene_list contains ENTREZ IDs as character strings
  gene_list <- as.character(na.omit(gene_list))

  if (is.null(gene_list) || length(gene_list) < 10) {
     log_message(paste("Skipping pathway analysis: Need at least 10 valid gene IDs, found", length(gene_list)))
    return(NULL)
  }

  OrgDb <- org.Hs.eg.db
  
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
if (!dir.exists(PARAMS$output_dir)) {
  dir.create(PARAMS$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", PARAMS$output_dir))
}

# Define output file paths
pca_plot_file <- file.path(PARAMS$output_dir, "pca_plot.pdf")
heatmap_file <- file.path(PARAMS$output_dir, "sample_correlation_heatmap.pdf")
top_degs_heatmap_file <- file.path(PARAMS$output_dir, "top_unique_degs_heatmap.pdf")
volcano_dir <- file.path(PARAMS$output_dir, "volcano_plots")
summary_dir <- file.path(PARAMS$output_dir, "summary_files")
comp_heatmap_dir <- file.path(PARAMS$output_dir, "comparison_heatmaps")

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
log_message(paste("Reading full metadata from:", PARAMS$metadata))
metadata_full <- read.csv(PARAMS$metadata, check.names = FALSE, row.names = "sample")

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
for (dds_file in list.files(PARAMS$dds_dir, pattern = "_vs_", include.dirs = TRUE, full.names = TRUE)) {
  comparison_dir <- basename(dds_file)
  if (dir.exists(dds_file) && file.exists(file.path(dds_file, "dds.rds"))) {
    dds_obj <- readRDS(file.path(dds_file, "dds.rds"))
    comparison_samples <- colnames(dds_obj)
    log_message(paste(" -", comparison_dir, ":", length(comparison_samples), "samples"))
    log_message(paste("   Samples:", paste(comparison_samples, collapse=", ")))
  }
}

# Read full count matrix
log_message(paste("Reading full count matrix from:", PARAMS$counts))
# Read the featureCounts output file
# Skip comment lines, use tab separator, specify header
raw_counts <- read.delim(PARAMS$counts, comment.char = "#", header = TRUE, sep = "\t", check.names = FALSE)

# Set Geneid as row names and select only count columns (7th column onwards)
# The first 6 columns are Geneid, Chr, Start, End, Strand, Length
rownames(raw_counts) <- raw_counts$Geneid
count_matrix_full <- raw_counts[, 7:ncol(raw_counts)]

# Clean up sample names in the count matrix (remove path prefixes if present)
# Example: results/star/C1/C1_Aligned.sortedByCoord.out.bam -> C1
colnames(count_matrix_full) <- sub(".*/(.*?)/\\1_Aligned.*", "\\1", colnames(count_matrix_full))
log_message(paste("Cleaned sample names:", paste(head(colnames(count_matrix_full)), collapse=", ")))


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

# Filter out low-count genes before VST (PARAMSional but recommended)
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ]

# Get variance stabilized transformed data for all samples
log_message("Performing variance stabilizing transformation on all samples")
vst_full <- vst(dds_full, blind = FALSE)

# Find all DESeq2 RDS files
log_message("Finding DESeq2 RDS files for comparisons")

# Search for comparison subdirectories within the specified dds_dir
comparisons <- list.dirs(PARAMS$dds_dir, full.names = TRUE, recursive = FALSE)
# Filter for directories that contain "_vs_" in their name
comparisons <- comparisons[grepl("_vs_", basename(comparisons))]

# Construct the full path to the dds.rds file within each comparison directory
dds_files <- file.path(comparisons, "dds.rds")

# Keep only the paths where the dds.rds file actually exists
dds_files <- dds_files[file.exists(dds_files)]

# Check if any files were found
if (length(dds_files) == 0) {
  stop("No DESeq2 RDS files found in comparison subdirectories within ", PARAMS$dds_dir, ". Comparison directories should contain '_vs_' in their name and a 'dds.rds' file.")
}

log_message(paste("Found", length(dds_files), "DESeq2 RDS file(s) to process"))
for (file in dds_files) {
  log_message(paste("  -", file))
}

# Create PCA plot using full vst data
log_message("Creating PCA plot using all samples")
pdf(pca_plot_file, width = PARAMS$pca_width, height = PARAMS$pca_height)

# The plotPCA function expects the intgroup to be present in colData(vst_full)
# Ensure that the 'condition' column is correctly placed in the colData
colData(vst_full)$condition <- metadata_full$condition

# Now generate the PCA data
pca_data <- plotPCA(vst_full, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = PARAMS$pca_point_size) +
  geom_text_repel(size = PARAMS$pca_text_size, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(PARAMS$pca_title) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Apply selected theme
theme_func <- match.fun(paste0("theme_", PARAMS$pca_theme))
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
num_colors_available <- brewer.pal.info[PARAMS$color_palette, "maxcolors"]
palette_colors <- brewer.pal(min(num_colors_needed, num_colors_available), PARAMS$color_palette)
if (num_colors_needed > num_colors_available) {
  warning(paste("Number of conditions exceeds colors in palette", PARAMS$color_palette, ". Colors will be recycled."))
  palette_colors <- rep(palette_colors, length.out = num_colors_needed)
}
condition_colors <- setNames(palette_colors, conditions_all)

# Create annotation data frame using full metadata
annotation_col <- data.frame(Condition = metadata_full$condition)
rownames(annotation_col) <- rownames(metadata_full)
ann_colors <- list(Condition = condition_colors)

pdf(heatmap_file, width = PARAMS$heatmap_width, height = PARAMS$heatmap_height)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = PARAMS$heatmap_title,
         fontsize = 10)
dev.off()
log_message(paste("Sample correlation heatmap saved to:", heatmap_file))

# Collect top DEGs from all comparisons
log_message("Collecting top DEGs from all comparisons")
all_top_degs <- list()
comparison_data <- list() # Store comparison DEGs and DDS path

for (dds_file in dds_files) {
  # Determine comparison name from the directory name
  comparison_name <- basename(dirname(dds_file))

  # Ensure comparison name looks valid (contains _vs_) - This check is now redundant due to list.dirs filtering
  # if (!grepl("_vs_", comparison_name)) {
  #     log_message(paste("Warning: Skipping file/directory", dds_file, "- does not appear to be a standard comparison name."))
  #     next # Skip to next dds_file
  # }

  log_message(paste("Processing DESeq2 object:", dds_file, "for comparison:", comparison_name))
  dds_obj <- readRDS(dds_file)

  # Get results
  res <- results(dds_obj)

  # Get top N DEGs by adjusted p-value for this comparison
  res_filtered <- as.data.frame(res) %>% filter(!is.na(padj))
  top_degs_comp <- rownames(res_filtered[order(res_filtered$padj), ])[1:min(PARAMS$top_n_degs, nrow(res_filtered))]

  all_top_degs <- c(all_top_degs, list(top_degs_comp))
  # Store DEGs and the correct path to the dds file
  comparison_data[[comparison_name]] <- list(degs = top_degs_comp, dds_path = dds_file)
}

# Get unique top DEGs across all comparisons
unique_top_degs <- unique(unlist(all_top_degs))
log_message(paste("Number of unique top DEGs across all comparisons:", length(unique_top_degs)))

# Create heatmap of unique top DEGs across all comparisons
log_message("Creating heatmap of unique top DEGs using only relevant samples")

for (comparison_name in names(comparison_data)) {
  top_degs_current_comp <- comparison_data[[comparison_name]]$degs # Get DEGs
  dds_file_path <- comparison_data[[comparison_name]]$dds_path     # Get correct DDS path

  if (length(top_degs_current_comp) > 0) {
    top_degs_current_comp_present <- intersect(top_degs_current_comp, rownames(assay(vst_full)))
    log_message(paste("Processing heatmap for:", comparison_name))
    log_message(paste("Found", length(top_degs_current_comp_present), "top DEGs for", comparison_name, "in VST data."))

    if (length(top_degs_current_comp_present) > 1) { # Need at least 2 genes to plot

      # Load the original DESeq2 object using the stored path
      log_message(paste("Loading DESeq2 object from:", dds_file_path))
      dds_comparison <- readRDS(dds_file_path) # Use the stored path
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
      heatmap_rownames <- get_gene_symbols(rownames(top_degs_counts_scaled), species = PARAMS$species)

      # Define file path for the heatmap
      # Ensure the subdirectory exists (it should have been created in the summary loop)
      comparison_heatmap_dir_specific <- file.path(comp_heatmap_dir, comparison_name)
      if (!dir.exists(comparison_heatmap_dir_specific)) { dir.create(comparison_heatmap_dir_specific, recursive = TRUE) }
      comp_heatmap_file_pdf <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_top_unique_degs_heatmap.pdf"))
      comp_heatmap_file_png <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_top_unique_degs_heatmap.png"))

      # Define color scale and breaks for the TOP DEGs heatmap (dynamic)
      # lim <- 2 # Set a fixed limit for the Z-score scale (e.g., -2 to +2) (Commented out)
      lim <- max(abs(top_degs_counts_scaled), na.rm = TRUE) # Use dynamic calculation based on top DEGs
      lim <- ceiling(lim * 10) / 10 # Round up slightly for better breaks
      if (lim == 0) lim <- 1 # Avoid zero range if data is flat
      my_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
      my_breaks <- seq(-lim, lim, length.out = 101) # Breaks based on dynamic lim

      # Generate heatmap using subsetted data and save as PDF
      pdf(comp_heatmap_file_pdf, width = PARAMS$top_degs_width, height = PARAMS$top_degs_height)
       pheatmap(top_degs_counts_scaled,
                  annotation_col = annotation_col_subset, # Use subsetted annotation
                  annotation_colors = ann_colors_subset,  # Use subset of annotation colors
                  labels_row = heatmap_rownames,
                  show_rownames = TRUE,
                  clustering_method = "ward.D2",
                  main = paste("Top Unique DEGs -", comparison_name), # Updated title
                  fontsize = 10,
                  fontsize_row = max(4, 10 - length(top_degs_current_comp_present) / 10),
                  color = my_colors,             # Apply custom color palette
                  breaks = my_breaks,            # Apply custom breaks
                  border_color = NA)             # Remove cell borders
      dev.off()
      log_message(paste("Heatmap of unique top DEGs for", comparison_name, "saved to:", comp_heatmap_file_pdf))

      # Generate heatmap using subsetted data and save as PNG
      png(comp_heatmap_file_png, width = PARAMS$top_degs_width, height = PARAMS$top_degs_height, units = "in", res = 300)
      pheatmap(top_degs_counts_scaled,
                  annotation_col = annotation_col_subset, # Use subsetted annotation
                  annotation_colors = ann_colors_subset,  # Use subset of annotation colors
                  labels_row = heatmap_rownames,
                  show_rownames = TRUE,
                  clustering_method = "ward.D2",
                  main = paste("Top Unique DEGs -", comparison_name), # Updated title
                  fontsize = 10,
                  fontsize_row = max(4, 10 - length(top_degs_current_comp_present) / 10),
                  color = my_colors,             # Apply custom color palette
                  breaks = my_breaks,            # Apply custom breaks
                  border_color = NA)             # Remove cell borders
      dev.off()
      log_message(paste("Heatmap of unique top DEGs for", comparison_name, "saved to:", comp_heatmap_file_png))

         # --- START: Add heatmap for all significant DEGs ---

         log_message(paste("Creating heatmap of ALL significant DEGs for:", comparison_name))

         # Get results again for this comparison to find all significant genes
         res_all <- results(dds_comparison)
         res_all_df <- as.data.frame(res_all) %>% filter(!is.na(padj))

         # Filter for significant genes (using volcano cutoffs for consistency in counts)
         sig_genes_all <- res_all_df %>% filter(padj < PARAMS$volcano_pval_cutoff)
         sig_gene_ids_all <- rownames(sig_genes_all)

         # Count up and down regulated among significant based on fold change cutoff
         up_regulated_count <- sum(sig_genes_all$log2FoldChange > PARAMS$volcano_fc_cutoff, na.rm = TRUE)
         down_regulated_count <- sum(sig_genes_all$log2FoldChange < -PARAMS$volcano_fc_cutoff, na.rm = TRUE)
         log_message(paste("Found", length(sig_gene_ids_all), "total significant genes (padj <", PARAMS$volcano_pval_cutoff, ").",
                           up_regulated_count, "up (LFC >", PARAMS$volcano_fc_cutoff,"),",
                           down_regulated_count, "down (LFC < -", PARAMS$volcano_fc_cutoff,")."))

         # Ensure these genes are present in the VST data
         sig_gene_ids_all_present <- intersect(sig_gene_ids_all, rownames(assay(vst_full)))
         log_message(paste("Found", length(sig_gene_ids_all_present), "significant DEGs present in VST data."))

         if (length(sig_gene_ids_all_present) > 1 && length(relevant_samples) > 1) {
             # Subset VST data for all significant genes and relevant samples
             all_sig_counts_subset <- assay(vst_full)[sig_gene_ids_all_present, relevant_samples, drop=FALSE]

             # Scale the subsetted counts
             all_sig_counts_scaled <- t(scale(t(all_sig_counts_subset)))
             all_sig_counts_scaled[is.nan(all_sig_counts_scaled)] <- 0 # Handle NaNs

             # Create title with unicode arrows
             heatmap_title_all_sig <- paste0("Heatmap of All Significant DEGs (p.adj < ", PARAMS$volcano_pval_cutoff, ")\n",
                                            comparison_name, "\n",
                                            up_regulated_count, " \u2191 and ", down_regulated_count, " \u2193",
                                            " (LFC cutoff: ", PARAMS$volcano_fc_cutoff, ")")

             # Define file paths
             all_sig_heatmap_file_pdf <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_all_significant_degs_heatmap.pdf"))
             all_sig_heatmap_file_png <- file.path(comparison_heatmap_dir_specific, paste0(comparison_name, "_all_significant_degs_heatmap.png"))

             # Define a FIXED color scale and breaks for the ALL significant DEGs heatmap
             lim_fixed <- 2 # Set the fixed limit
             if (lim_fixed == 0) lim_fixed <- 1 # Avoid zero range
             # Reuse the same color palette function if desired, or define a new one
             my_colors_fixed <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
             my_breaks_fixed <- seq(-lim_fixed, lim_fixed, length.out = 101) # Breaks based on fixed lim

             # Generate heatmap (PDF) - No row names as there are too many
             # Adjust height/width if needed, maybe make height larger for many genes? Using same for now.
             pdf(all_sig_heatmap_file_pdf, width = PARAMS$top_degs_width, height = PARAMS$top_degs_height)
              pheatmap(all_sig_counts_scaled,
                      annotation_col = annotation_col_subset,
                      annotation_colors = ann_colors_subset,
                      color = my_colors_fixed,   # Use fixed colors
                      breaks = my_breaks_fixed,  # Use fixed breaks
                      border_color = NA,
                      main = heatmap_title_all_sig,
                      show_rownames = FALSE, # Hide row names for clarity
                      show_colnames = TRUE,
                      clustering_method = "ward.D2",
                      fontsize = 10,
                      fontsize_title = 12) # Slightly larger title font
             dev.off()
             log_message(paste("Heatmap of all significant DEGs for", comparison_name, "saved to:", all_sig_heatmap_file_pdf))

             # Generate heatmap (PNG) - No row names
             png(all_sig_heatmap_file_png, width = PARAMS$top_degs_width, height = PARAMS$top_degs_height, units = "in", res = 300)
              pheatmap(all_sig_counts_scaled,
                      annotation_col = annotation_col_subset,
                      annotation_colors = ann_colors_subset,
                      color = my_colors_fixed,   # Use fixed colors
                      breaks = my_breaks_fixed,  # Use fixed breaks
                      border_color = NA,
                      main = heatmap_title_all_sig,
                      show_rownames = FALSE, # Hide row names for clarity
                      show_colnames = TRUE,
                      clustering_method = "ward.D2",
                      fontsize = 10,
                      fontsize_title = 12) # Slightly larger title font
             dev.off()
             log_message(paste("Heatmap of all significant DEGs for", comparison_name, "saved to:", all_sig_heatmap_file_png))

         } else {
             log_message(paste("Skipping heatmap of all significant DEGs for", comparison_name, ": Not enough significant genes (", length(sig_gene_ids_all_present), ") or samples (", length(relevant_samples), ") found."))
         }

         # --- END: Add heatmap for all significant DEGs ---

       # Removed the redundant 'else' block that logged "DESeq2 object not found"
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
   res_df$gene_symbol <- get_gene_symbols(res_df$gene_id, species = PARAMS$species)

   # Calculate summary statistics
   total_genes <- nrow(res_df)
   significant_genes <- sum(res_df$padj < PARAMS$volcano_pval_cutoff & !is.na(res_df$padj))
   up_regulated <- sum(res_df$padj < PARAMS$volcano_pval_cutoff & res_df$log2FoldChange > PARAMS$volcano_fc_cutoff & !is.na(res_df$padj))
   down_regulated <- sum(res_df$padj < PARAMS$volcano_pval_cutoff & res_df$log2FoldChange < -PARAMS$volcano_fc_cutoff & !is.na(res_df$padj))

   # Get top up and down regulated genes
   sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < PARAMS$volcano_pval_cutoff, ]

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
   sig_genes_df <- res_df[res_df$padj < PARAMS$volcano_pval_cutoff & !is.na(res_df$padj), ]
   sig_genes_df <- sig_genes_df[order(sig_genes_df$padj), ]
   write.csv(sig_genes_df, sig_genes_file, row.names = FALSE)
   log_message(paste("Significant genes saved to:", sig_genes_file))

   # Create up and down regulated gene lists
   up_genes_file <- file.path(comparison_summary_dir, "up_regulated.csv")
   up_genes_df <- res_df[res_df$padj < PARAMS$volcano_pval_cutoff &
                         res_df$log2FoldChange > PARAMS$volcano_fc_cutoff &
                         !is.na(res_df$padj), ]
   up_genes_df <- up_genes_df[order(up_genes_df$log2FoldChange, decreasing = TRUE), ]
   write.csv(up_genes_df, up_genes_file, row.names = FALSE)
   log_message(paste("Up-regulated genes saved to:", up_genes_file))

   down_genes_file <- file.path(comparison_summary_dir, "down_regulated.csv")
   down_genes_df <- res_df[res_df$padj < PARAMS$volcano_pval_cutoff &
                           res_df$log2FoldChange < -PARAMS$volcano_fc_cutoff &
                           !is.na(res_df$padj), ]
   down_genes_df <- down_genes_df[order(down_genes_df$log2FoldChange), ]
   write.csv(down_genes_df, down_genes_file, row.names = FALSE)
   log_message(paste("Down-regulated genes saved to:", down_genes_file))

   # Run pathway analysis if requested
   if (PARAMS$run_pathway) {
     log_message("Running pathway analysis")

     # Get ENTREZ IDs for significant genes
     sig_entrez <- get_entrez_ids(sig_genes_df$gene_id, species = PARAMS$species)
     sig_entrez <- sig_entrez[!is.na(sig_entrez)]

     # Run pathway analysis
     pathway_results <- run_pathway_analysis(sig_entrez, species = PARAMS$species)

     if (!is.null(pathway_results) && nrow(pathway_results) > 0) {
       pathway_file <- file.path(comparison_summary_dir, "pathway_analysis.csv")
       write.csv(pathway_results, pathway_file, row.names = FALSE)
       log_message(paste("Pathway analysis results saved to:", pathway_file))
     } else {
       log_message("No significant pathway enrichment found")
     }

     # Separate pathway analysis for up and down regulated genes
     up_entrez <- get_entrez_ids(up_genes_df$gene_id, species = PARAMS$species)
     up_entrez <- up_entrez[!is.na(up_entrez)]

     if (length(up_entrez) >= 10) {
       up_pathway_results <- run_pathway_analysis(up_entrez, species = PARAMS$species)

       if (!is.null(up_pathway_results) && nrow(up_pathway_results) > 0) {
         up_pathway_file <- file.path(comparison_summary_dir, "up_pathway_analysis.csv")
         write.csv(up_pathway_results, up_pathway_file, row.names = FALSE)
         log_message(paste("Up-regulated pathway analysis results saved to:", up_pathway_file))
       }
     }

     down_entrez <- get_entrez_ids(down_genes_df$gene_id, species = PARAMS$species)
     down_entrez <- down_entrez[!is.na(down_entrez)]

     if (length(down_entrez) >= 10) {
       down_pathway_results <- run_pathway_analysis(down_entrez, species = PARAMS$species)

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
   res_df$gene_symbol <- get_gene_symbols(res_df$gene_id, species = PARAMS$species)

   # Add significance column
   res_df$significant <- "Not Significant"
   res_df$significant[res_df$padj < PARAMS$volcano_pval_cutoff & !is.na(res_df$padj) & res_df$log2FoldChange > PARAMS$volcano_fc_cutoff] <- "Up-regulated"
   res_df$significant[res_df$padj < PARAMS$volcano_pval_cutoff & !is.na(res_df$padj) & res_df$log2FoldChange < -PARAMS$volcano_fc_cutoff] <- "Down-regulated"
   res_df$significant <- factor(res_df$significant, levels = c("Up-regulated", "Down-regulated", "Not Significant"))

   # Define colors for volcano plot
   volcano_colors <- c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")

   # Get top genes for labeling (based on p-value and fold change)
   top_genes <- res_df %>%
     filter(!is.na(padj) & padj < PARAMS$volcano_pval_cutoff & abs(log2FoldChange) > PARAMS$volcano_fc_cutoff) %>%
     mutate(rank_metric = -log10(padj) * abs(log2FoldChange)) %>%
     arrange(desc(rank_metric)) %>%
     head(PARAMS$volcano_top_n)

   # Create volcano plot
   volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
     geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
     geom_hline(yintercept = -log10(PARAMS$volcano_pval_cutoff), linetype = "dashed") +
     geom_vline(xintercept = c(-PARAMS$volcano_fc_cutoff, PARAMS$volcano_fc_cutoff), linetype = "dashed") +
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
     pdf(volcano_file, width = PARAMS$volcano_width, height = PARAMS$volcano_height)
     print(volcano_plot)
     dev.off()
     log_message(paste("Volcano plot saved to:", volcano_file))
   }, error = function(e) {
     warning(paste("Failed to save PDF volcano plot:", e$message))
   })

   # Also save as PNG for easier viewing - with error handling
   volcano_png <- file.path(comparison_volcano_dir, "volcano_plot.png")
   tryCatch({
     ggsave(volcano_png, volcano_plot, width = PARAMS$volcano_width, height = PARAMS$volcano_height, dpi = 300)
     log_message(paste("Volcano plot saved to:", volcano_png))
   }, error = function(e) {
     warning(paste("Failed to save PNG volcano plot:", e$message))
   })
 }

 log_message("Visualization generation complete")
