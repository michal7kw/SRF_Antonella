#!/bin/bash
# shellcheck disable=SC2086,SC2046

# ==============================================================================
# Script: Generate Heatmaps around TSS for CUT&Tag Data
#
# Description:
#   Calculates a signal matrix and plots a heatmap of CUT&Tag signal
#   (from BigWig files) centered around gene Transcription Start Sites (TSSs).
#   Uses deepTools computeMatrix and plotHeatmap.
#   Processes one sample at a time based on a command-line index.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Define input paths (BigWig, Gene Annotation).
#   3. Calculate signal matrix around TSS using computeMatrix.
#   4. Generate heatmap plot using plotHeatmap.
#   5. Clean up temporary files.
#
# Usage:
#   bash 6_heatmap_tss.sh <sample_index>
# Example:
#   bash 6_heatmap_tss.sh 0  # Processes the first sample (GFP_1)
#
# Parallel Execution:
#   Use GNU Parallel or xargs for parallel processing across samples:
#   N=6 # Number of parallel jobs
#   seq 0 5 | parallel -j $N bash 6_heatmap_tss.sh {}
#   # OR
#   # printf "%s\n" {0..5} | xargs -I {} -P $N bash 6_heatmap_tss.sh {}
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
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/6_heatmap_tss"
LOG_DIR="${SCRIPT_DIR}/logs/6_heatmap_tss"
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_heatmap_tss" # Base for temporary files

# Sample names array (should match BigWig filenames without .bw)
SAMPLES=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# deepTools computeMatrix Parameters
COMPUTE_MATRIX_UPSTREAM=5000      # Base pairs upstream of TSS
COMPUTE_MATRIX_DOWNSTREAM=5000    # Base pairs downstream of TSS
COMPUTE_MATRIX_BINSIZE=50         # Bin size for averaging signal
COMPUTE_MATRIX_THREADS=4          # Number of threads for computeMatrix

# deepTools plotHeatmap Parameters (Adjust as needed)
PLOT_HEATMAP_COLORMAP="Blues"        # Colormap (e.g., Blues, viridis, RdYlBu)
PLOT_HEATMAP_SORT_REGIONS="descend"  # Sort regions (no, descend, ascend, mean, median, min, max, sum)
PLOT_HEATMAP_SORT_USING="mean"       # Sort using metric (mean, median, max, min, sum, std, region_length)
PLOT_HEATMAP_ZMIN="0"                # Minimum value for color scale (optional, remove for auto)
PLOT_HEATMAP_ZMAX="10"               # Maximum value for color scale (optional, remove for auto)
PLOT_HEATMAP_HEIGHT=20               # Heatmap height in cm
PLOT_HEATMAP_WIDTH=8                 # Heatmap width in cm
PLOT_HEATMAP_XAXIS_LABEL="Distance from TSS (bp)" # X-axis label
PLOT_HEATMAP_REGIONS_LABEL="Genes" # Label for the regions (rows)
# PLOT_HEATMAP_SAMPLES_LABEL will be set per sample below
# PLOT_HEATMAP_TITLE will be set per sample below

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
    check_tool plotHeatmap # Added check for plotHeatmap
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
    local heatmap_plot_png="${output_dir}/${sample_name}_tss_heatmap.png" # Renamed variable
    local heatmap_plot_pdf="${output_dir}/${sample_name}_tss_heatmap.pdf" # Renamed variable

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

    # Step 2: Plot Heatmap from Matrix
    log "Step 2: Running plotHeatmap for ${sample_name}..."
    local plotHeatmap_log="${log_dir}/plotHeatmap.log"
    local PLOT_HEATMAP_TITLE="Signal Heatmap around TSS (${sample_name})" # Define title here

    # Plot PNG
    plotHeatmap \
        -m "${matrix_gz}" \
        -out "${heatmap_plot_png}" \
        --colorMap ${PLOT_HEATMAP_COLORMAP} \
        --sortRegions ${PLOT_HEATMAP_SORT_REGIONS} \
        --sortUsing ${PLOT_HEATMAP_SORT_USING} \
        --whatToShow "heatmap and colorbar" \
        --zMin ${PLOT_HEATMAP_ZMIN} \
        --zMax ${PLOT_HEATMAP_ZMAX} \
        --heatmapHeight ${PLOT_HEATMAP_HEIGHT} \
        --heatmapWidth ${PLOT_HEATMAP_WIDTH} \
        --xAxisLabel "${PLOT_HEATMAP_XAXIS_LABEL}" \
        --regionsLabel "${PLOT_HEATMAP_REGIONS_LABEL}" \
        --samplesLabel "${sample_name}" \
        --plotTitle "${PLOT_HEATMAP_TITLE}" \
        --verbose \
        >& "${plotHeatmap_log}" \
        || die "plotHeatmap (PNG) failed for ${sample_name}. Check log: ${plotHeatmap_log}"

    # Plot PDF (reuse log, append)
    plotHeatmap \
        -m "${matrix_gz}" \
        -out "${heatmap_plot_pdf}" \
        --colorMap ${PLOT_HEATMAP_COLORMAP} \
        --sortRegions ${PLOT_HEATMAP_SORT_REGIONS} \
        --sortUsing ${PLOT_HEATMAP_SORT_USING} \
        --whatToShow "heatmap and colorbar" \
        --zMin ${PLOT_HEATMAP_ZMIN} \
        --zMax ${PLOT_HEATMAP_ZMAX} \
        --heatmapHeight ${PLOT_HEATMAP_HEIGHT} \
        --heatmapWidth ${PLOT_HEATMAP_WIDTH} \
        --xAxisLabel "${PLOT_HEATMAP_XAXIS_LABEL}" \
        --regionsLabel "${PLOT_HEATMAP_REGIONS_LABEL}" \
        --samplesLabel "${sample_name}" \
        --plotTitle "${PLOT_HEATMAP_TITLE}" \
        --verbose \
        >> "${plotHeatmap_log}" 2>&1 \
        || die "plotHeatmap (PDF) failed for ${sample_name}. Check log: ${plotHeatmap_log}"


    # Verify plotHeatmap output
    [[ -s "$heatmap_plot_png" ]] || die "plotHeatmap failed to produce PNG output file: $heatmap_plot_png. Check log: ${plotHeatmap_log}"
    [[ -s "$heatmap_plot_pdf" ]] || die "plotHeatmap failed to produce PDF output file: $heatmap_plot_pdf. Check log: ${plotHeatmap_log}"
    log "plotHeatmap complete for ${sample_name}."

    # --- Final Summary ---
    log "--- Summary for ${sample_name} ---"
    log "Output Directory: ${output_dir}"
    log "Matrix File: ${matrix_gz}"
    log "Heatmap Plot (PNG): ${heatmap_plot_png}"
    log "Heatmap Plot (PDF): ${heatmap_plot_pdf}"
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