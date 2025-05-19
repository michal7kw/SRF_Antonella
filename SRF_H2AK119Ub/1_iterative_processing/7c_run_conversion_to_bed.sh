#!/bin/bash

# Documentation:
# This script converts the differential_peaks.csv file (output from DiffBind analysis)
# into a BED file suitable for visualization in the Integrative Genomics Viewer (IGV).

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
#   "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_lfc_significant.csv"
#   "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_significant.csv"
#   "analysis/7_differential_binding_strict_v2/differential_peaks.csv"
)

OUTPUT_BED_FILENAMES=(
  "analysis/7_differential_binding_v2/differential_peaks_fdr_lfc_significant.bed"
  "analysis/7_differential_binding_v2/differential_peaks_fdr_significant.bed"
  "analysis/7_differential_binding_v2/differential_peaks.bed"
#   "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_lfc_significant.bed"
#   "analysis/7_differential_binding_strict_v2/differential_peaks_fdr_significant.bed"
#   "analysis/7_differential_binding_strict_v2/differential_peaks.bed"
)

# R script for conversion
R_SCRIPT_PATH="scripts/7c_convert_diffpeaks_to_bed.R"

# Parameters for the R script
# Choose the column for the BED score. Common choices:
# "Fold" (log2 fold change)
# "FDR" (False Discovery Rate)
# "p.value" (p-value)
SCORE_COLUMN="Fold" # Defaulting to Fold change

# Set to "true" if your IGV genome uses UCSC-style chromosome names (e.g., "chr1", "chrM")
# Set to "false" if your IGV genome uses non-UCSC style (e.g., "1", "MT")
# The 7_differential_binding_v2.R script standardizes to non-UCSC style (e.g., "1", "X").
# So, if your IGV uses hg38 from UCSC, you might need to add "chr" prefix.
ADD_CHR_PREFIX="true" # Ensure UCSC style output to match GTF

log_message "Starting batch conversion of differential peaks CSV files to BED format for IGV."

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    log_message "ERROR: Rscript command not found. Please ensure R is installed and in your PATH."
    exit 1
fi

# Check if the R conversion script exists
if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    log_message "ERROR: R conversion script not found at $R_SCRIPT_PATH"
    exit 1
fi

# Loop through the defined input/output pairs
# Ensure INPUT_CSV_FILENAMES is not empty before proceeding
if [ ${#INPUT_CSV_FILENAMES[@]} -eq 0 ]; then
    log_message "ERROR: No input CSV files defined in INPUT_CSV_FILENAMES array."
    exit 1
fi
if [ ${#INPUT_CSV_FILENAMES[@]} -ne ${#OUTPUT_BED_FILENAMES[@]} ]; then
    log_message "ERROR: Mismatch between the number of input CSV files and output BED files."
    exit 1
fi

for i in "${!INPUT_CSV_FILENAMES[@]}"; do
    INPUT_CSV_PATH="${INPUT_CSV_FILENAMES[$i]}"
    OUTPUT_BED_PATH="${OUTPUT_BED_FILENAMES[$i]}"

    log_message "---------------------------------------------------------------------"
    log_message "Input CSV: ${INPUT_CSV_PATH}"
    log_message "Output BED: ${OUTPUT_BED_PATH}"

    # Check if the input CSV file exists
    if [[ ! -f "$INPUT_CSV_PATH" ]]; then
        log_message "WARNING: Input CSV file not found at $INPUT_CSV_PATH. Skipping this pair."
        continue # Skip to the next iteration
    fi

    log_message "R script: ${R_SCRIPT_PATH}"
    log_message "Score column for BED: ${SCORE_COLUMN}"
    log_message "Add 'chr' prefix: ${ADD_CHR_PREFIX}"

    # Construct Rscript command
    R_COMMAND="Rscript ${R_SCRIPT_PATH} \
        --input_csv \"${INPUT_CSV_PATH}\" \
        --output_bed \"${OUTPUT_BED_PATH}\" \
        --score_column \"${SCORE_COLUMN}\""

    if [[ "$ADD_CHR_PREFIX" == "true" ]]; then
        R_COMMAND="${R_COMMAND} --add_chr_prefix"
    fi

    # Execute the R script
    log_message "Running R script for conversion for ${INPUT_CSV_PATH}..."
    if eval "$R_COMMAND"; then
        log_message "Successfully converted ${INPUT_CSV_PATH} to ${OUTPUT_BED_PATH}"
    else
        log_message "ERROR: R script execution failed for ${INPUT_CSV_PATH}. Check output for details."
        # Decide if you want to exit on first error or continue with other files
        # exit 1 # Uncomment to exit on first error
    fi
done

log_message "---------------------------------------------------------------------"
log_message "Batch conversion script finished."