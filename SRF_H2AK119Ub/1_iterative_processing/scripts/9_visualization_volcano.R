# This script creates a volcano plot for differential peaks analysis results.
# It helps interpret differential binding analysis, focusing on promoter-associated peaks.
#
# Input:
#   - A CSV file containing promoter-associated gene annotations.
#     Required columns: 'fold_change', 'FDR', 'SYMBOL'.
#
# Output:
#   - A PDF file named 'volcano_plot_promoters.pdf' in a subdirectory:
#     {output_dir}/volcano_plot_promoters.pdf
#
# Command-line arguments:
#   1. output_dir: Base directory for output files.
#   2. input_csv_file: Path to the input CSV file with promoter gene data.
#   3. peak_type: A character string (e.g., "broad") used in naming the output path.
#
# Dependencies:
#   - R packages: ggplot2, ggrepel
#   - A 'scripts/utils.R' file providing log_message() and save_plot() functions.

# Load necessary libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
})

# Source helper functions (ensure this path is correct relative to script execution)
source("scripts/utils.R")

#' Create and save a volcano plot for promoter peak analysis.
#'
#' @param output_dir Base directory where the plot's subdirectory will be created.
#' @param promoter_data A data frame with promoter-associated gene annotations.
#'                      Must contain 'fold_change', 'FDR', and 'SYMBOL' columns.
#' @return The ggplot object for the volcano plot.
create_and_save_volcano_plot <- function(output_dir, promoter_data, file_name) {

    # Validate input data
    required_cols <- c("fold_change", "FDR", "SYMBOL")
    if (!all(required_cols %in% names(promoter_data))) {
        missing_cols <- required_cols[!required_cols %in% names(promoter_data)]
        stop(sprintf("Input data frame 'promoter_data' is missing required columns: %s",
                     paste(missing_cols, collapse = ", ")))
    }

    # Prepare output directory
    plot_subdir <- file.path(output_dir)
    create_dirs(plot_subdir) # Assumes create_dirs is from utils.R and handles multiple paths

    # Handle FDR values of 0 to prevent -log10(0) = Inf issues
    if (any(promoter_data$FDR == 0, na.rm = TRUE)) {
        log_message("Found FDR values of 0. Replacing with a very small number for -log10 transformation.")
        min_positive_fdr <- min(promoter_data$FDR[promoter_data$FDR > 0 & is.finite(promoter_data$FDR)], na.rm = TRUE)
        replacement_fdr <- if (is.finite(min_positive_fdr) && min_positive_fdr > 0) {
            min_positive_fdr * 0.01 
        } else {
            .Machine$double.eps # A very small positive number
        }
        promoter_data$FDR[promoter_data$FDR == 0 & !is.na(promoter_data$FDR)] <- replacement_fdr
    }
    
    # Calculate significant gene counts for plot annotation
    promoter_data$is_significant <- promoter_data$FDR < 0.01 & !is.na(promoter_data$FDR)
    total_sig_count <- sum(promoter_data$is_significant, na.rm = TRUE)
    sig_yaf_count <- sum(promoter_data$fold_change > 0 & promoter_data$is_significant, na.rm = TRUE)
    sig_gfp_count <- sum(promoter_data$fold_change < 0 & promoter_data$is_significant, na.rm = TRUE)

    # Prepare data for labeling top genes
    # Filter out NA FDRs before ordering for top genes
    promoter_data_valid_fdr <- promoter_data[!is.na(promoter_data$FDR) & is.finite(promoter_data$FDR), ]
    
    top_yaf_genes <- head(promoter_data_valid_fdr[promoter_data_valid_fdr$fold_change > 0, ][order(promoter_data_valid_fdr[promoter_data_valid_fdr$fold_change > 0, "FDR"]), ], 8)
    top_gfp_genes <- head(promoter_data_valid_fdr[promoter_data_valid_fdr$fold_change < 0, ][order(promoter_data_valid_fdr[promoter_data_valid_fdr$fold_change < 0, "FDR"]), ], 8)
    top_sig_genes_to_label <- rbind(top_yaf_genes, top_gfp_genes)

    # Determine y-axis minimum for annotation text placement
    y_values_for_plot <- -log10(promoter_data$FDR)
    y_min_finite <- min(y_values_for_plot[is.finite(y_values_for_plot)], na.rm = TRUE)
    annotation_y_pos <- if (is.finite(y_min_finite)) y_min_finite + 0.5 else 0.5


    # Create the volcano plot
    volcano_plot <- ggplot(promoter_data, aes(x = fold_change, y = -log10(FDR))) +
        geom_point(aes(color = is_significant, size = abs(fold_change)), alpha = 0.6, na.rm = TRUE) +
        scale_color_manual(
            values = c("TRUE" = "red3", "FALSE" = "grey50"),
            name = "FDR < 0.01",
            labels = c("TRUE" = "Yes", "FALSE" = "No"),
            drop = FALSE # Ensure all levels are shown even if data is filtered
        ) +
        scale_size_continuous(
            range = c(0.5, 2.5),
            limits = c(0, max(abs(promoter_data$fold_change), na.rm = TRUE)),
            name = "abs(Fold Change)"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey50", alpha = 0.5) +
        annotate("text",
                 x = 0,
                 y = annotation_y_pos,
                 label = sprintf("Significant (FDR < 0.01):\nTotal: %d | YAF-enriched: %d | GFP-enriched: %d",
                               total_sig_count, sig_yaf_count, sig_gfp_count),
                 hjust = 0.5, vjust = 0, size = 3.0, color = "black", fontface = "bold", lineHeight = 0.9) +
        geom_text_repel(
            data = top_sig_genes_to_label,
            aes(label = SYMBOL),
            size = 2.5,
            box.padding = 0.5,
            point.padding = 0.2,
            force = 8,
            max.overlaps = 15,
            segment.color = "grey50",
            segment.alpha = 0.6,
            na.rm = TRUE
        ) +
        labs(
            x = "log2 Fold Change (YAF/GFP)",
            y = "-log10(FDR)",
            title = file_name
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5),
            axis.title = element_text(size = rel(1.1)),
            legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey90"),
            plot.margin = margin(10, 10, 10, 10, unit = "pt")
        )

    # Save the plot
    output_file_path <- file.path(plot_subdir, paste0(file_name, "_volcano_plot.pdf"))
    save_plot(volcano_plot, output_file_path, width = 9, height = 8) # Assumes save_plot is from utils.R
    log_message(sprintf("Volcano plot saved to: %s", output_file_path))

    return(volcano_plot)
}

# --- Main script execution ---
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 3) {
        stop("Usage: Rscript 9_visualization_volcano.R <output_dir> <input_csv_file> <file_name>", call. = FALSE)
    }

    output_dir_arg <- args[1]
    input_csv_file_arg <- args[2]
    file_name <- args[3]

    log_message(sprintf("Starting volcano plot generation"))
    log_message(sprintf("Output directory: %s", output_dir_arg))
    log_message(sprintf("Input CSV: %s", input_csv_file_arg))

    if (!file.exists(input_csv_file_arg)) {
        stop(sprintf("Input file not found: %s", input_csv_file_arg))
    }

    promoter_data_df <- read.csv(input_csv_file_arg, stringsAsFactors = FALSE)

    # Create and save the volcano plot
    create_and_save_volcano_plot(output_dir_arg, promoter_data_df, file_name)

    log_message("Volcano plot generation process completed successfully.")
}

# Run the main function when script is executed
if (sys.nframe() == 0) {
    main()
}