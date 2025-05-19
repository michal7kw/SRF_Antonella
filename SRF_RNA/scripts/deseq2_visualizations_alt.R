#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
})

# Set working directory to the specified path

results_path <- "/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP"
output_path <- "/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/viz_deseq2_alt/YAF_vs_GFP"


PARAMS <- list(
  dds_dir = results_path,
  metadata = "/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_RNA/metadata.csv",
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
  pca_point_size = 4,
  pca_text_size = 4,
  pca_theme = "bw",
  color_palette = "Set1",
  heatmap_color = "RdBu",
  show_gene_names = TRUE,
  max_gene_names = 30,
  clustering_method = "ward.D2",
  pca_components = "1,2",
  verbose = TRUE
)

# Function to print messages if verbose is TRUE
log_message <- function(message) {
  if (PARAMS$verbose) {
    cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  }
}

log_message("Starting enhanced visualization generation")

# Create output directory if it doesn't exist
if (!dir.exists(PARAMS$output_dir)) {
  dir.create(PARAMS$output_dir, recursive = TRUE)
  log_message(paste("Created output directory:", PARAMS$output_dir))
}

# Define output file paths
pca_plot_file <- file.path(PARAMS$output_dir, "pca_plot.pdf")
heatmap_file <- file.path(PARAMS$output_dir, "sample_correlation_heatmap.pdf")
top_degs_heatmap_file <- file.path(PARAMS$output_dir, "top_degs_heatmap.pdf")

# Read metadata
log_message(paste("Reading metadata from:", PARAMS$metadata))
metadata <- read.csv(PARAMS$metadata, check.names = FALSE)

# Find all DESeq2 RDS files
log_message(paste0("Looking for DESeq2 RDS file (dds.rds) directly in: ", PARAMS$dds_dir))
potential_dds_file <- file.path(PARAMS$dds_dir, "dds.rds")

if (file.exists(potential_dds_file)) {
  dds_files <- c(potential_dds_file) # Store as a character vector, suitable for the existing loops
  log_message(paste0("Found dds.rds at: ", potential_dds_file))
} else {
  dds_files <- character(0) # Initialize as empty if not found, to trigger the error message in the next block
  log_message(paste0("dds.rds NOT FOUND directly in: ", PARAMS$dds_dir, ". Please ensure the file 'dds.rds' exists at this specific location."))
}

if (length(dds_files) == 0) {
  stop("No DESeq2 RDS files found in ", PARAMS$dds_dir)
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

# Parse PCA components
pca_comps <- as.numeric(strsplit(PARAMS$pca_components, ",")[[1]])
if (length(pca_comps) != 2) {
  stop("PCA components must be specified as two comma-separated numbers")
}

# Create PCA plot
log_message("Creating PCA plot")
pdf(pca_plot_file, width = PARAMS$pca_width, height = PARAMS$pca_height)

# Get PCA data for all components
pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create enhanced PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = PARAMS$pca_point_size) +
  geom_text_repel(size = PARAMS$pca_text_size, show.legend = FALSE, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(PARAMS$pca_title) +
  theme_get() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Apply selected theme
if (PARAMS$pca_theme == "bw") {
  pca_plot <- pca_plot + theme_bw()
} else if (PARAMS$pca_theme == "classic") {
  pca_plot <- pca_plot + theme_classic()
} else if (PARAMS$pca_theme == "minimal") {
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
condition_colors <- brewer.pal(min(length(unique(metadata$condition)), 9), PARAMS$color_palette)
names(condition_colors) <- unique(metadata$condition)
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- metadata$sample
ann_colors <- list(Condition = condition_colors)

# Choose color palette for heatmap
if (PARAMS$heatmap_color == "viridis") {
  color_func <- viridis(100)
} else if (PARAMS$heatmap_color == "magma") {
  color_func <- viridis(100, PARAMSion = "A")
} else if (PARAMS$heatmap_color == "inferno") {
  color_func <- viridis(100, PARAMSion = "B")
} else if (PARAMS$heatmap_color == "plasma") {
  color_func <- viridis(100, PARAMSion = "C")
} else {
  color_func <- colorRampPalette(rev(brewer.pal(9, PARAMS$heatmap_color)))(100)
}

pdf(heatmap_file, width = PARAMS$heatmap_width, height = PARAMS$heatmap_height)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = PARAMS$heatmap_title,
         fontsize = 10,
         color = color_func,
         clustering_method = PARAMS$clustering_method)
dev.off()
log_message(paste("Sample correlation heatmap saved to:", heatmap_file))

# Collect top DEGs from all comparisons
log_message("Collecting top DEGs from all comparisons")
all_top_degs <- list()

for (dds_file in dds_files) {
  comparison_name <- basename(dirname(dds_file))
  log_message(paste("Processing DESeq2 object:", comparison_name))
  dds_obj <- readRDS(dds_file)
  
  # Get results
  res <- results(dds_obj)
  
  # Store results for later use
  
  # Get top N DEGs by adjusted p-value
  top_degs <- rownames(res[order(res$padj), ])[1:PARAMS$top_n_degs]
  all_top_degs <- c(all_top_degs, list(top_degs))
}

# Get unique top DEGs
unique_top_degs <- unique(unlist(all_top_degs))
log_message(paste("Number of unique top DEGs:", length(unique_top_degs)))

# Create heatmap of top DEGs
log_message("Creating heatmap of top DEGs")

if (length(unique_top_degs) == 0) {
  log_message("No unique top DEGs found. Skipping heatmap generation for top DEGs.")
} else {
  # Extract normalized counts for top DEGs
  # Ensure unique_top_degs are present in vst assay rownames before subsetting
  valid_top_degs <- unique_top_degs[unique_top_degs %in% rownames(assay(vst))]
  if (length(valid_top_degs) == 0) {
    log_message("None of the unique top DEGs were found in the VST assay. Skipping heatmap generation.")
  } else {
    if (length(valid_top_degs) < length(unique_top_degs)) {
      log_message(paste("Warning:", length(unique_top_degs) - length(valid_top_degs),
                        "DEGs were not found in the VST assay and will be excluded from the heatmap."))
    }
    top_degs_counts <- assay(vst)[valid_top_degs, , drop = FALSE] # drop=FALSE ensures it remains a matrix

    # Scale the counts for better visualization
    # Check for constant rows before scaling to avoid errors with scale()
    # A row is constant if all its values are the same (e.g., all zeros if a gene has no counts in any sample)
    # Such rows would result in NaNs after scaling if sd is 0.
    row_sds <- apply(top_degs_counts, 1, sd)
    non_constant_rows <- row_sds > 1e-6 # Check for standard deviation greater than a small epsilon
    
    if(sum(non_constant_rows) == 0) {
        log_message("All top DEGs have constant expression across samples (e.g. all zeros). Cannot scale data for heatmap. Skipping heatmap.")
    } else {
        if(sum(!non_constant_rows) > 0) {
            log_message(paste("Warning:", sum(!non_constant_rows),
                              "DEGs have constant expression across samples and will be excluded from scaling and the heatmap if they are the only ones."))
        }
        # Proceed with scaling only for non-constant rows
        top_degs_counts_to_scale <- top_degs_counts[non_constant_rows, , drop = FALSE]
        
        if(nrow(top_degs_counts_to_scale) > 0) {
            top_degs_counts_scaled <- t(scale(t(top_degs_counts_to_scale)))

            # Determine whether to show row names
            show_rownames_heatmap <- PARAMS$show_gene_names
            if (PARAMS$show_gene_names && nrow(top_degs_counts_scaled) > PARAMS$max_gene_names) {
              log_message(paste("Number of DEGs (", nrow(top_degs_counts_scaled), ") exceeds max_gene_names (", PARAMS$max_gene_names, "). Hiding gene names in heatmap.", sep=""))
              show_rownames_heatmap <- FALSE
            }

            # Create the heatmap
            pdf(top_degs_heatmap_file, width = PARAMS$top_degs_width, height = PARAMS$top_degs_height)
            pheatmap(top_degs_counts_scaled, # Use scaled data
                     annotation_col = annotation_col,
                     annotation_colors = ann_colors,
                     show_rownames = show_rownames_heatmap, # Use defined variable
                     clustering_method = PARAMS$clustering_method,
                     main = PARAMS$top_degs_title, # Use parameter for title
                     fontsize = 10,
                     fontsize_row = 6,
                     color = color_func)
            dev.off()
            log_message(paste("Top DEGs heatmap saved to:", top_degs_heatmap_file))
        } else {
            log_message("No non-constant DEGs to plot after filtering. Skipping heatmap generation.")
        }
    }
  }
}

log_message("Enhanced visualization generation complete")

