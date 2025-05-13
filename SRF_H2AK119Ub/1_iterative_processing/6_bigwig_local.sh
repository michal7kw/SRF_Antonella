#!/bin/bash
# shellcheck disable=SC2086,SC2046

# ==============================================================================
# Script: BigWig Generation and Library Complexity Estimation (Local)
#
# Description:
#   Generates a normalized bigWig coverage file and estimates library complexity
#   from a deduplicated BAM file for a single sample.
#   This script is adapted for local execution, removing SLURM dependencies.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Check for required tools (bamCoverage, picard/java).
#   3. Validate input BAM file existence.
#   4. Generate normalized bigWig file using deepTools bamCoverage.
#   5. Estimate library complexity using Picard EstimateLibraryComplexity.
#   6. Validate output files.
#   7. Clean up temporary files.
#
# Usage:
#   bash 6_bigwig_local.sh <sample_index>
# Example:
#   bash 6_bigwig_local.sh 0  # Processes the first sample (GFP_1)
#
# Parallel Execution:
#   Use GNU Parallel or xargs for parallel processing across samples:
#   N=6 # Number of parallel jobs (adjust based on local resources)
#   seq 0 5 | parallel -j $N bash 6_bigwig_local.sh {}
#   # OR
#   # printf "%s\n" {0..5} | xargs -I {} -P $N bash 6_bigwig_local.sh {}
#
# Requirements:
#   - Conda environment with deepTools and Picard Tools installed and activated.
#   - Input BAM files (deduplicated) and BAI index files.
# ==============================================================================

# --- Configuration ---
# Exit on error, treat unset variables as errors, propagate exit status through pipes
set -euo pipefail

# Script directory for relative paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- Parameters & Paths (Modify as needed) ---
# Input data directory structure (relative to SCRIPT_DIR)
ALIGNMENT_DIR="${SCRIPT_DIR}/analysis/3_alignment"
# Output directory (relative to SCRIPT_DIR)
OUTPUT_DIR="${SCRIPT_DIR}/analysis/6_bigwig"
LOG_DIR="${SCRIPT_DIR}/logs/6_bigwig"
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_bigwig" # Base for temporary files

# Sample names array
SAMPLES=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# Number of processors for bamCoverage (adjust based on local machine)
NUM_PROCESSORS=8

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
    check_tool bamCoverage
    check_tool picard # Or check_tool java if using java -jar picard.jar
    check_tool mkdir
    check_tool date
    check_tool rm
    check_tool printf # Used for logging/output

    # --- Setup Directories & Paths ---
    local log_dir="${LOG_DIR}/${sample_name}" # Sample-specific log dir
    local timestamp=$(date +%Y%m%d%H%M%S%N) # Added nanoseconds for higher uniqueness
    # Define SAMPLE_TMP_DIR globally for trap cleanup
    # Ensure TMP_BASE_DIR exists before creating sample temp dir
    mkdir -p "${TMP_BASE_DIR}" || die "Failed to create base temporary directory: ${TMP_BASE_DIR}"
    SAMPLE_TMP_DIR="${TMP_BASE_DIR}/${sample_name}_${timestamp}"

    mkdir -p "${OUTPUT_DIR}" || die "Failed to create output directory: ${OUTPUT_DIR}"
    mkdir -p "${log_dir}" || die "Failed to create log directory: ${log_dir}"
    mkdir -p "${SAMPLE_TMP_DIR}" || die "Failed to create sample temporary directory: ${SAMPLE_TMP_DIR}"
    # Set TMPDIR for tool internal use if needed
    export TMPDIR="${SAMPLE_TMP_DIR}"

    # --- Trap for Cleanup ---
    # Ensure SAMPLE_TMP_DIR is cleaned up on exit, error, interrupt, or termination
    # Needs SAMPLE_TMP_DIR to be defined *before* trap is set.
    trap cleanup EXIT ERR INT TERM

    # --- Define File Paths ---
    local input_bam="${ALIGNMENT_DIR}/${sample_name}.dedup.bam"
    local input_bai="${ALIGNMENT_DIR}/${sample_name}.dedup.bam.bai"
    local output_bw="${OUTPUT_DIR}/${sample_name}.bw"
    local output_complexity="${OUTPUT_DIR}/${sample_name}_complexity.txt"
    local bamcov_log="${log_dir}/bamCoverage.log"
    local picard_log="${log_dir}/picard_complexity.log"

    # --- Verify Input Files ---
    [[ -f "$input_bam" ]] || die "Input BAM not found: $input_bam"
    [[ -f "$input_bai" ]] || die "Input BAI not found: $input_bai"

    # --- Activate Conda Environment ---
    # Adjust conda path and environment name if necessary
    # log "Activating conda environment..."
    # source /opt/common/tools/ric.cosr/miniconda3/bin/activate # Adjust path if needed
    # conda activate snakemake # Adjust environment name if needed
    # log "Conda environment activated."

    # --- Pipeline Steps ---

    # Step 1: Generate BigWig file
    log "Step 1: Generating bigWig file for ${sample_name}..."
    bamCoverage \
        -b "$input_bam" \
        -o "$output_bw" \
        --binSize 10 \
        --normalizeUsing RPKM \
        --smoothLength 30 \
        --extendReads \
        --centerReads \
        --numberOfProcessors ${NUM_PROCESSORS} \
        >& "$bamcov_log" \
        || die "bamCoverage failed for ${sample_name}. Check log: $bamcov_log"

    # Verify bamCoverage output
    [[ -s "$output_bw" ]] || die "bamCoverage failed to produce output file: $output_bw. Check log: $bamcov_log"
    log "BigWig generation complete for ${sample_name}."

    # Step 2: Estimate Library Complexity
    log "Step 2: Calculating library complexity for ${sample_name}..."
    picard EstimateLibraryComplexity \
        I="$input_bam" \
        O="$output_complexity" \
        TMP_DIR="$SAMPLE_TMP_DIR" \
        >& "$picard_log" \
        || die "Picard EstimateLibraryComplexity failed for ${sample_name}. Check log: $picard_log"

    # Verify Picard output
    [[ -s "$output_complexity" ]] || die "Picard failed to produce output file: $output_complexity. Check log: $picard_log"
    log "Library complexity estimation complete for ${sample_name}."

    # --- Final Summary ---
    log "--- Summary for ${sample_name} ---"
    log "Output BigWig: ${output_bw}"
    log "Output Complexity: ${output_complexity}"
    log "Log Directory: ${log_dir}"
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