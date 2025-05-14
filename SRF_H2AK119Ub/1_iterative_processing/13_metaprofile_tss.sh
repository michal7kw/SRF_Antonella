#!/bin/bash
# shellcheck disable=SC2086,SC2046

# ==============================================================================
# Script: Generate Metaprofiles around TSS for CUT&Tag Data
#
# Description:
#   Calculates and plots the average signal profile (metaprofile) of
#   CUT&Tag signal (from BigWig files) centered around gene Transcription
#   Start Sites (TSSs). Uses deepTools computeMatrix and plotProfile.
#   Processes one sample at a time based on a command-line index.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Define input paths (BigWig, Gene Annotation).
#   3. Calculate signal matrix around TSS using computeMatrix.
#   4. Generate metaprofile plot using plotProfile.
#   5. Clean up temporary files.
#
# Usage:
#   bash 6_metaprofile_tss.sh <sample_index>
# Example:
#   bash 6_metaprofile_tss.sh 0  # Processes the first sample (GFP_1)
#
# Parallel Execution:
#   Use GNU Parallel or xargs for parallel processing across samples:
#   N=6 # Number of parallel jobs
#   seq 0 5 | parallel -j $N bash 6_metaprofile_tss.sh {}
#   # OR
#   # printf "%s\n" {0..5} | xargs -I {} -P $N bash 6_metaprofile_tss.sh {}
#
# Requirements:
#   - Conda environment with deepTools installed and activated.
#   - Input BigWig files.
#   - Gene annotation file (GTF or BED format).
# ==============================================================================

# --- Configuration ---
# Exit on error, treat unset variables as errors, propagate exit status through pipes
set -euo pipefail

# Script directory for relative paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- Parameters & Paths (Modify as needed) ---

# === USER INPUT REQUIRED ===
# Path to the gene annotation file (GTF or BED format)
# Example: GENE_ANNOTATION_GTF="${SCRIPT_DIR}/../../COMMON_DATA/hg38.refGene.gtf"
GENE_ANNOTATION_GTF="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"

# Directory containing the input BigWig files. Assumes files are named {sample_name}.bw
# Example: BIGWIG_DIR="${SCRIPT_DIR}/analysis/4_bigwig"
BIGWIG_DIR="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/6_bigwig"
# === END USER INPUT ===

# Output directory (relative to SCRIPT_DIR)
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/13_metaprofile_tss"
LOG_DIR="${SCRIPT_DIR}/logs/13_metaprofile_tss"
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_metaprofile_tss" # Base for temporary files

# Sample names array (should match BigWig filenames without .bw)
SAMPLES=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# deepTools Parameters
COMPUTE_MATRIX_UPSTREAM=5000      # Base pairs upstream of TSS
COMPUTE_MATRIX_DOWNSTREAM=5000    # Base pairs downstream of TSS
COMPUTE_MATRIX_BINSIZE=50         # Bin size for averaging signal
COMPUTE_MATRIX_THREADS=4          # Number of threads for computeMatrix
PLOT_PROFILE_COLORS="blue"        # Color for the profile line (can be space-separated list for multiple samples/groups)
PLOT_PROFILE_HEIGHT=15            # Plot height in cm
PLOT_PROFILE_WIDTH=8              # Plot width in cm

# --- Helper Functions ---
log() {
    # Logs messages to stderr with a timestamp
    printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*" >&2
}

die() {
    # Logs an error message to stderr and exits
    log "ERROR: $*"
    exit 1
}

check_tool() {
    # Checks if a command exists in the PATH
    command -v "$1" >/dev/null 2>&1 || die "Required tool '$1' not found in PATH. Ensure conda environment is active and contains the tool."
}

cleanup() {
    # Removes the temporary directory if it exists
    # Declared globally so trap can access it
    log "Executing cleanup..."
    if [[ -n "${SAMPLE_TMP_DIR:-}" && -d "${SAMPLE_TMP_DIR}" ]]; then
        log "Removing temporary directory: ${SAMPLE_TMP_DIR}"
        rm -rf "${SAMPLE_TMP_DIR}"
        log "Temporary directory removed."
    else
         log "Skipping cleanup: Temporary directory path invalid ('${SAMPLE_TMP_DIR:-}') or does not exist."
    fi
}

# --- Main Script Logic ---

main() {
    # --- Argument Parsing & Validation ---
    if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
        printf "Usage: %s <sample_index>\n" "$(basename "$0")" >&2
        printf "Provide the integer index (0-%d) of the sample to process.\n" $((${#SAMPLES[@]} - 1)) >&2
        exit 1
    fi
    local sample_index=$1

    if [[ $sample_index -lt 0 ]] || [[ $sample_index -ge ${#SAMPLES[@]} ]]; then
        die "Invalid sample index $sample_index. Must be between 0 and $((${#SAMPLES[@]} - 1))."
    fi
    local sample_name=${SAMPLES[$sample_index]}
    log "Starting processing for sample: ${sample_name} (Index: ${sample_index})"

    # --- Check Dependencies ---
    check_tool computeMatrix
    check_tool plotProfile
    check_tool mkdir
    check_tool date
    check_tool rm
    check_tool printf # Used for logging/output

    # --- Setup Directories & Paths ---
    local output_dir="${BASE_OUTPUT_DIR}/${sample_name}" # Sample-specific output dir
    local log_dir="${LOG_DIR}/${sample_name}" # Sample-specific log dir
    local timestamp=$(date +%Y%m%d%H%M%S%N) # Added nanoseconds for higher uniqueness
    # Define SAMPLE_TMP_DIR globally for trap cleanup
    # Ensure TMP_BASE_DIR exists before creating sample temp dir
    mkdir -p "${TMP_BASE_DIR}" || die "Failed to create base temporary directory: ${TMP_BASE_DIR}"
    SAMPLE_TMP_DIR="${TMP_BASE_DIR}/${sample_name}_${timestamp}"

    mkdir -p "${output_dir}" || die "Failed to create output directory: ${output_dir}"
    mkdir -p "${log_dir}" || die "Failed to create log directory: ${log_dir}"
    mkdir -p "${SAMPLE_TMP_DIR}" || die "Failed to create sample temporary directory: ${SAMPLE_TMP_DIR}"
    # Set TMPDIR for potential tool use, though not strictly needed by deepTools usually
    export TMPDIR="${SAMPLE_TMP_DIR}"

    # --- Trap for Cleanup ---
    # Ensure SAMPLE_TMP_DIR is cleaned up on exit, error, interrupt, or termination
    # Needs SAMPLE_TMP_DIR to be defined *before* trap is set.
    trap cleanup EXIT ERR INT TERM

    # --- Define File Paths ---
    local input_bw="${BIGWIG_DIR}/${sample_name}.bw"
    local matrix_gz="${output_dir}/${sample_name}_tss_matrix.mat.gz"
    local matrix_tab="${output_dir}/${sample_name}_tss_matrix.tab" # Optional matrix values
    local matrix_bed="${output_dir}/${sample_name}_tss_matrix.bed" # Optional region info
    local profile_plot_png="${output_dir}/${sample_name}_tss_profile.png"
    local profile_plot_pdf="${output_dir}/${sample_name}_tss_profile.pdf" # Also generate PDF

    # --- Verify Input Files/Placeholders ---
    [[ -f "$input_bw" ]] || die "Input BigWig not found: $input_bw"
    if [[ "$GENE_ANNOTATION_GTF" == "PLEASE_PROVIDE_PATH/hg38.refGene.gtf" ]]; then
        die "Gene annotation file path (GENE_ANNOTATION_GTF) is a placeholder. Please update the script with the correct path."
    fi
    [[ -f "$GENE_ANNOTATION_GTF" ]] || die "Gene annotation file not found: $GENE_ANNOTATION_GTF"
    if [[ "$BIGWIG_DIR" == "PLEASE_PROVIDE_PATH/analysis/4_bigwig" ]]; then
        die "BigWig directory path (BIGWIG_DIR) is a placeholder. Please update the script with the correct path."
    fi
    [[ -d "$BIGWIG_DIR" ]] || die "BigWig directory not found: $BIGWIG_DIR"

    # --- Pipeline Steps ---

    # Step 1: Compute Matrix around TSS
    log "Step 1: Running computeMatrix reference-point for ${sample_name}..."
    local computeMatrix_log="${log_dir}/computeMatrix.log"
    computeMatrix reference-point \
        --referencePoint TSS \
        -b ${COMPUTE_MATRIX_UPSTREAM} \
        -a ${COMPUTE_MATRIX_DOWNSTREAM} \
        -R "${GENE_ANNOTATION_GTF}" \
        -S "${input_bw}" \
        -o "${matrix_gz}" \
        --outFileNameMatrix "${matrix_tab}" \
        --outFileSortedRegions "${matrix_bed}" \
        --skipZeros \
        --missingDataAsZero \
        --binSize ${COMPUTE_MATRIX_BINSIZE} \
        --numberOfProcessors ${COMPUTE_MATRIX_THREADS} \
        --verbose \
        >& "${computeMatrix_log}" \
        || die "computeMatrix failed for ${sample_name}. Check log: ${computeMatrix_log}"

    # Verify computeMatrix output
    [[ -s "$matrix_gz" ]] || die "computeMatrix failed to produce output file: $matrix_gz. Check log: ${computeMatrix_log}"
    log "computeMatrix complete for ${sample_name}."

    # Step 2: Plot Profile from Matrix
    log "Step 2: Running plotProfile for ${sample_name}..."
    local plotProfile_log="${log_dir}/plotProfile.log"

    # Plot PNG
    plotProfile \
        -m "${matrix_gz}" \
        --plotType lines \
        -out "${profile_plot_png}" \
        --samplesLabel "${sample_name}" \
        --colors ${PLOT_PROFILE_COLORS} \
        --plotHeight ${PLOT_PROFILE_HEIGHT} \
        --plotWidth ${PLOT_PROFILE_WIDTH} \
        --yAxisLabel "Average Signal" \
        --plotTitle "Signal Profile around TSS (${sample_name})" \
        --verbose \
        >& "${plotProfile_log}" \
        || die "plotProfile (PNG) failed for ${sample_name}. Check log: ${plotProfile_log}"

    # Plot PDF (reuse log, append)
    plotProfile \
        -m "${matrix_gz}" \
        --plotType lines \
        -out "${profile_plot_pdf}" \
        --samplesLabel "${sample_name}" \
        --colors ${PLOT_PROFILE_COLORS} \
        --plotHeight ${PLOT_PROFILE_HEIGHT} \
        --plotWidth ${PLOT_PROFILE_WIDTH} \
        --yAxisLabel "Average Signal" \
        --plotTitle "Signal Profile around TSS (${sample_name})" \
        --verbose \
        >> "${plotProfile_log}" 2>&1 \
        || die "plotProfile (PDF) failed for ${sample_name}. Check log: ${plotProfile_log}"


    # Verify plotProfile output
    [[ -s "$profile_plot_png" ]] || die "plotProfile failed to produce PNG output file: $profile_plot_png. Check log: ${plotProfile_log}"
    [[ -s "$profile_plot_pdf" ]] || die "plotProfile failed to produce PDF output file: $profile_plot_pdf. Check log: ${plotProfile_log}"
    log "plotProfile complete for ${sample_name}."

    # --- Final Summary ---
    log "--- Summary for ${sample_name} ---"
    log "Output Directory: ${output_dir}"
    log "Matrix File: ${matrix_gz}"
    log "Profile Plot (PNG): ${profile_plot_png}"
    log "Profile Plot (PDF): ${profile_plot_pdf}"
    log "Temporary files were in: ${SAMPLE_TMP_DIR}" # Will be cleaned up by trap
    log "--- Processing Complete for ${sample_name} ---"

} # End of main function

# --- Execute Main ---
# Make SAMPLE_TMP_DIR global and initialize for the trap
SAMPLE_TMP_DIR=""
# Pass all script arguments to the main function
main "$@"

# Cleanup is handled by the trap on EXIT/ERR/INT/TERM
log "Script finished successfully."
exit 0

# --- Example Run Commands ---
# Single sample:
# bash 6_metaprofile_tss.sh 0

# Parallel using GNU Parallel:
# N=6; seq 0 5 | parallel -j $N bash 6_metaprofile_tss.sh {}

# Parallel using xargs:
# N=6; printf "%s\n" {0..5} | xargs -I {} -P $N bash 6_metaprofile_tss.sh {}