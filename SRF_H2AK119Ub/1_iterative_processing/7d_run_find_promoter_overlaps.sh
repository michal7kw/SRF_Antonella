#!/bin/bash

# Documentation:
# This script identifies differential peaks that overlap with gene promoter regions.
# Promoter regions are defined based on a GTF annotation file.

# Error handling
set -e  # Exit immediately if a command exits with non-zero status
set -u  # Treat unset variables as an error
set -o pipefail  # Fail pipeline if any command fails

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Define working directory (assuming script is run from 1_iterative_processing)
WORKDIR="."
cd "$WORKDIR" || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define input and output file arrays
INPUT_CSV_FILENAMES=(
    "analysis/7_differential_binding_v2/differential_peaks_fdr_lfc_significant.csv"
    "analysis/7_differential_binding_v2/differential_peaks_fdr_significant.csv"
    "analysis/7_differential_binding_v2/differential_peaks.csv"
    # "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_lfc_significant.csv"
    # "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_significant.csv"
    # "analysis/7_differential_binding_strict_v2/differential_peaks.csv"
)

# Gene annotation GTF file
GTF_FILE_PATH="../../COMMON_DATA/gencode.v43.basic.annotation.gtf"

# R script for conversion
R_SCRIPT_PATH="scripts/7d_find_promoter_overlapping_diffpeaks.R"

# Parameters for the R script
PROMOTER_UPSTREAM=2000  # Distance upstream of TSS (bp)
PROMOTER_DOWNSTREAM=500   # Distance downstream of TSS (bp)
SCORE_COLUMN="Fold"       # Column from CSV for BED score (e.g., Fold, FDR)

# Chromosome naming style for inputs and output.
# Your 7_differential_binding_v2.R script standardizes peak files to non-UCSC (e.g., "1", "X").
# Gencode GTF files typically use UCSC style ("chr1", "chrX").
# The R script 9_... will attempt to reconcile these.
# Set --add_chr_prefix_output if you want "chr" in the final BED file for IGV.
PEAK_CHR_STYLE="NCBI" # e.g. 1, 2, X, MT (as produced by your step 7 script)
GTF_CHR_STYLE="UCSC"  # e.g. chr1, chr2, chrX, chrM (common for Gencode)
ADD_CHR_PREFIX_OUTPUT="true" # Ensure UCSC style output to match GTF

log_message "Starting batch identification of differential peaks in promoter regions."

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    log_message "ERROR: Rscript command not found. Please ensure R is installed and in your PATH."
    exit 1
fi

# Check if the R script exists
if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    log_message "ERROR: R script not found at $R_SCRIPT_PATH"
    exit 1
fi

# Check if the GTF file exists
if [[ ! -f "$GTF_FILE_PATH" ]]; then
    log_message "ERROR: GTF annotation file not found at $GTF_FILE_PATH"
    log_message "Please verify the path to 'gencode.v43.basic.annotation.gtf'."
    exit 1
fi
log_message "Using GTF Annotation: ${GTF_FILE_PATH}"
log_message "Using R script: ${R_SCRIPT_PATH}"
log_message "Promoter definition: ${PROMOTER_UPSTREAM}bp upstream, ${PROMOTER_DOWNSTREAM}bp downstream of TSS"
log_message "Score column for BED: ${SCORE_COLUMN}"
log_message "Peak chromosome style (input): ${PEAK_CHR_STYLE}"
log_message "GTF chromosome style (input): ${GTF_CHR_STYLE}"
log_message "Add 'chr' prefix to output BED: ${ADD_CHR_PREFIX_OUTPUT}"

# Ensure INPUT_CSV_FILENAMES is not empty before proceeding
if [ ${#INPUT_CSV_FILENAMES[@]} -eq 0 ]; then
    log_message "ERROR: No input CSV files defined in INPUT_CSV_FILENAMES array."
    exit 1
fi

# Loop through the defined input CSV files
for CURRENT_INPUT_CSV_FILENAME in "${INPUT_CSV_FILENAMES[@]}"; do
    # Derive output BED filename
    BASE_NAME=$(basename "${CURRENT_INPUT_CSV_FILENAME}" .csv)
    CURRENT_OUTPUT_BED_FILENAME="${BASE_NAME}_in_promoters.bed"

    # Define paths inside the loop
    INPUT_CSV_PATH="${CURRENT_INPUT_CSV_FILENAME}"
    OUTPUT_BED_PATH="${CURRENT_OUTPUT_BED_FILENAME}"

    log_message "---------------------------------------------------------------------"
    log_message "Processing Input CSV: ${CURRENT_INPUT_CSV_FILENAME}"

    # Check if the input CSV file exists
    if [[ ! -f "$INPUT_CSV_PATH" ]]; then
        log_message "WARNING: Input CSV file not found at $INPUT_CSV_PATH. Skipping this file."
        continue # Skip to the next iteration
    fi

    log_message "Output BED will be: ${OUTPUT_BED_PATH}" # Clarify output path before R script call

    # Construct Rscript command
    R_COMMAND="Rscript ${R_SCRIPT_PATH} \
        --diff_peaks_csv \"${INPUT_CSV_PATH}\" \
        --gtf_file \"${GTF_FILE_PATH}\" \
        --output_bed \"${OUTPUT_BED_PATH}\" \
        --promoter_upstream ${PROMOTER_UPSTREAM} \
        --promoter_downstream ${PROMOTER_DOWNSTREAM} \
        --score_column \"${SCORE_COLUMN}\" \
        --peak_chr_style \"${PEAK_CHR_STYLE}\" \
        --gtf_chr_style \"${GTF_CHR_STYLE}\""

    if [[ "$ADD_CHR_PREFIX_OUTPUT" == "true" ]]; then
        R_COMMAND="${R_COMMAND} --add_chr_prefix_output"
    fi

    # Execute the R script
    log_message "Running R script to find peaks in promoters for ${CURRENT_INPUT_CSV_FILENAME}..."
    if eval "$R_COMMAND"; then
        log_message "Successfully found promoter overlaps for ${CURRENT_INPUT_CSV_FILENAME}. Output: ${OUTPUT_BED_PATH}"
    else
        log_message "ERROR: R script execution failed for ${CURRENT_INPUT_CSV_FILENAME}. Check output for details."
        # Decide if you want to exit on first error or continue
        # exit 1 # Uncomment to exit on first error
    fi
done

log_message "---------------------------------------------------------------------"
log_message "Batch promoter overlap analysis finished."