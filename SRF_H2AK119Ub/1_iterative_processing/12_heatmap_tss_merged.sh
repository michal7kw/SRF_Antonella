#!/bin/bash
# shellcheck disable=SC2086,SC2046,SC2207

# ==============================================================================
# Script: Generate Merged Heatmaps around TSS for CUT&Tag Data per Sample Group
#
# Description:
#   Calculates a signal matrix and plots a heatmap of CUT&Tag signal
#   (from multiple BigWig files corresponding to replicas of a sample group)
#   centered around gene Transcription Start Sites (TSSs).
#   The signals from replicas are averaged by computeMatrix.
#   Uses deepTools computeMatrix and plotHeatmap.
#   Processes one sample group at a time based on a command-line index.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Define input paths (BigWig directory, Gene Annotation).
#   3. Identify replica BigWig files for the specified sample group.
#   4. Calculate a single signal matrix around TSS using computeMatrix with all replica BigWigs.
#   5. Generate a single heatmap plot using plotHeatmap.
#   6. Clean up temporary files.
#
# Usage:
#   bash 13_heatmap_tss_merged.sh <sample_group_index>
# Example:
#   bash 13_heatmap_tss_merged.sh 0  # Processes the first sample group (e.g., GFP)
#
# Parallel Execution:
#   Define SAMPLE_GROUP_PREFIXES array first. If N_GROUPS is the number of elements:
#   N_JOBS=$((${#SAMPLE_GROUP_PREFIXES[@]})) # Number of parallel jobs, one per group
#   seq 0 $((${N_JOBS} - 1)) | parallel -j ${N_JOBS} bash 13_heatmap_tss_merged.sh {}
#
# Requirements:
#   - Conda environment with deepTools installed and activated.
#   - Input BigWig files, named like {sample_group_prefix}_{replica_id}.bw (e.g., GFP_1.bw).
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

# Directory containing the input BigWig files.
# Assumes files are named {sample_group_prefix}_{replica_id}.bw (e.g., GFP_1.bw, YAF_1.bw)
# Example: BIGWIG_DIR="${SCRIPT_DIR}/analysis/10_bigwig"
BIGWIG_DIR="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/10_bigwig"
# === END USER INPUT ===

# Output directory (relative to SCRIPT_DIR) - Adjusted for merged analysis
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/12_heatmap_tss/merged"
LOG_DIR="${SCRIPT_DIR}/logs/12_heatmap_tss_merged"
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_heatmap_tss_merged" # Base for temporary files

# Sample group prefixes (should match the part of BigWig filenames before replica ID)
SAMPLE_GROUP_PREFIXES=(GFP YAF) # e.g., GFP will group GFP_1.bw, GFP_2.bw, etc.

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
# PLOT_HEATMAP_SAMPLES_LABEL will be set per sample group below
# PLOT_HEATMAP_TITLE will be set per sample group below

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
    if [[ -n "${SAMPLE_GROUP_TMP_DIR:-}" && -d "${SAMPLE_GROUP_TMP_DIR}" ]]; then
        log "Removing temporary directory: ${SAMPLE_GROUP_TMP_DIR}"
        rm -rf "${SAMPLE_GROUP_TMP_DIR}"
        log "Temporary directory removed."
    else
         log "Skipping cleanup: Temporary directory path invalid ('${SAMPLE_GROUP_TMP_DIR:-}') or does not exist."
    fi
}

# --- Main Script Logic ---

main() {
    # --- Argument Parsing & Validation ---
    if [[ $# -ne 1 || ! "$1" =~ ^[0-9]+$ ]]; then
        printf "Usage: %s <sample_group_index>\n" "$(basename "$0")" >&2
        printf "Provide the integer index (0-%d) of the sample group to process from SAMPLE_GROUP_PREFIXES.\n" $((${#SAMPLE_GROUP_PREFIXES[@]} - 1)) >&2
        exit 1
    fi
    local sample_group_index=$1

    if [[ $sample_group_index -lt 0 ]] || [[ $sample_group_index -ge ${#SAMPLE_GROUP_PREFIXES[@]} ]]; then
        die "Invalid sample group index $sample_group_index. Must be between 0 and $((${#SAMPLE_GROUP_PREFIXES[@]} - 1))."
    fi
    local sample_group_prefix=${SAMPLE_GROUP_PREFIXES[$sample_group_index]}
    
    # Define output_dir and final plot paths early for the check
    local output_dir="${BASE_OUTPUT_DIR}/${sample_group_prefix}"
    local final_heatmap_plot_png="${output_dir}/${sample_group_prefix}_tss_heatmap.png"

    # Check if final plot outputs already exist
    if [[ -f "${final_heatmap_plot_png}" && -s "${final_heatmap_plot_png}" ]]; then
        log "Final heatmap PNG already exists for ${sample_group_prefix}:"
        log "  - ${final_heatmap_plot_png}"
        log "Skipping processing for this group."
        exit 0
    fi
    
    log "Starting processing for sample group: ${sample_group_prefix} (Index: ${sample_group_index})"
    # --- Check Dependencies ---
    check_tool computeMatrix
    check_tool plotHeatmap
    check_tool bigwigCompare # Added for averaging BigWigs
    check_tool mkdir
    check_tool date
    check_tool rm
    check_tool printf
    check_tool find # Used to find replica BigWig files
    check_tool sort # Used to sort replica BigWig files

    # --- Setup Directories & Paths ---
    local output_dir="${BASE_OUTPUT_DIR}/${sample_group_prefix}" # Sample-group-specific output dir
    local log_file_dir="${LOG_DIR}/${sample_group_prefix}" # Sample-group-specific log dir
    local timestamp=$(date +%Y%m%d%H%M%S%N) # Added nanoseconds for higher uniqueness
    # Define SAMPLE_GROUP_TMP_DIR globally for trap cleanup
    # Ensure TMP_BASE_DIR exists before creating sample group temp dir
    mkdir -p "${TMP_BASE_DIR}" || die "Failed to create base temporary directory: ${TMP_BASE_DIR}"
    SAMPLE_GROUP_TMP_DIR="${TMP_BASE_DIR}/${sample_group_prefix}_${timestamp}"

    mkdir -p "${output_dir}" || die "Failed to create output directory: ${output_dir}"
    mkdir -p "${log_file_dir}" || die "Failed to create log directory: ${log_file_dir}"
    mkdir -p "${SAMPLE_GROUP_TMP_DIR}" || die "Failed to create sample group temporary directory: ${SAMPLE_GROUP_TMP_DIR}"
    # Set TMPDIR for potential tool use, though not strictly needed by deepTools usually
    export TMPDIR="${SAMPLE_GROUP_TMP_DIR}"

    # --- Trap for Cleanup ---
    # Ensure SAMPLE_GROUP_TMP_DIR is cleaned up on exit, error, interrupt, or termination
    trap cleanup EXIT ERR INT TERM

    # --- Define File Paths ---
    # Find all BigWig files for the current sample group prefix
    local input_bw_files=()
    # Use find and sort to get a consistent list of replica files
    # The glob pattern ensures we only pick files like GFP_1.bw, GFP_anything.bw
    local found_files
    found_files=$(find "${BIGWIG_DIR}" -maxdepth 1 -name "${sample_group_prefix}_*.bw" -print0 | sort -zV | xargs -0)

    if [[ -z "$found_files" ]]; then
        die "No BigWig files found for prefix '${sample_group_prefix}' in directory '${BIGWIG_DIR}' with pattern '${sample_group_prefix}_*.bw'"
    fi
    # Convert space-separated string from xargs back to array
    read -r -a input_bw_files <<< "$found_files"

    log "Found ${#input_bw_files[@]} BigWig files for sample group ${sample_group_prefix}:"
    for bwf in "${input_bw_files[@]}"; do
        log "  - ${bwf}"
        [[ -f "$bwf" ]] || die "Input BigWig file not found: $bwf (listed for group ${sample_group_prefix})"
    done

    local matrix_gz="${output_dir}/${sample_group_prefix}_tss_matrix.mat.gz"
    local matrix_tab="${output_dir}/${sample_group_prefix}_tss_matrix.tab" # Optional matrix values
    local matrix_bed="${output_dir}/${sample_group_prefix}_tss_matrix.bed" # Optional region info
    local heatmap_plot_png="${output_dir}/${sample_group_prefix}_tss_heatmap.png"

    # --- Verify Input Files/Placeholders ---
    # Individual BigWig file checks are done above after finding them.
    if [[ "$GENE_ANNOTATION_GTF" == "PLEASE_PROVIDE_PATH/hg38.refGene.gtf" ]]; then # Keep this check as an example
        die "Gene annotation file path (GENE_ANNOTATION_GTF) is a placeholder. Please update the script with the correct path."
    fi
    [[ -f "$GENE_ANNOTATION_GTF" ]] || die "Gene annotation file not found: $GENE_ANNOTATION_GTF"
    if [[ "$BIGWIG_DIR" == "PLEASE_PROVIDE_PATH/analysis/4_bigwig" ]]; then # Keep this check as an example
        die "BigWig directory path (BIGWIG_DIR) is a placeholder. Please update the script with the correct path."
    fi
    [[ -d "$BIGWIG_DIR" ]] || die "BigWig directory not found: $BIGWIG_DIR"

    # --- Pipeline Steps ---

    local num_replicas=${#input_bw_files[@]}
    local effective_bw_for_matrix # This will be the single BigWig input for computeMatrix
    local matrix_scale_factor=1.0 # Scale factor for computeMatrix
    local final_summed_bw_path="${SAMPLE_GROUP_TMP_DIR}/${sample_group_prefix}_averaged_sum.bw" # Standardized name

    if [[ $num_replicas -eq 0 ]]; then
        die "Logic error: No BigWig files found for ${sample_group_prefix} before processing."
    elif [[ $num_replicas -eq 1 ]]; then
        log "Single replica found for ${sample_group_prefix}. Using it directly."
        effective_bw_for_matrix="${input_bw_files[0]}"
    else # num_replicas > 1
        log "Multiple replicas ($num_replicas) found for ${sample_group_prefix}."
        if [[ -f "${final_summed_bw_path}" && -s "${final_summed_bw_path}" ]]; then
            log "Found existing summed BigWig: ${final_summed_bw_path}. Using it."
            effective_bw_for_matrix="${final_summed_bw_path}"
        else
            log "No existing summed BigWig found or it is empty. Summing replicas using bigwigCompare..."
            local bigwigCompare_sum_log="${log_file_dir}/bigwigCompare_sum.log"
            >"${bigwigCompare_sum_log}" # Ensure log is fresh

            local sum_in_progress_bw="${input_bw_files[0]}"
            local temp_sum_storage="${SAMPLE_GROUP_TMP_DIR}/sum_intermediates"
            mkdir -p "${temp_sum_storage}"

            for (( i=1; i<num_replicas; i++ )); do
                local file_to_add="${input_bw_files[$i]}"
                local current_step_output_bw
                if [[ $i -eq $((num_replicas - 1)) ]]; then # Last iteration, save to final path
                    current_step_output_bw="${final_summed_bw_path}"
                else
                    current_step_output_bw="${temp_sum_storage}/${sample_group_prefix}_sum_intermediate_step_${i}.bw"
                fi
                
                log "Summing (step $((i)) of $((num_replicas - 1))): ${sum_in_progress_bw} + ${file_to_add} -> ${current_step_output_bw}"
                bigwigCompare \
                    --bigwig1 "${sum_in_progress_bw}" \
                    --bigwig2 "${file_to_add}" \
                    --operation add \
                    --binSize 10 \
                    --skipNonCoveredRegions \
                    --numberOfProcessors "${COMPUTE_MATRIX_THREADS}" \
                    -o "${current_step_output_bw}" \
                    >> "${bigwigCompare_sum_log}" 2>&1 \
                    || die "bigwigCompare add operation failed. Check log: ${bigwigCompare_sum_log}"
                [[ -s "${current_step_output_bw}" ]] || die "bigwigCompare add operation did not produce output file: ${current_step_output_bw}"
                
                sum_in_progress_bw="${current_step_output_bw}"
            done
            effective_bw_for_matrix="${final_summed_bw_path}"
            log "All $num_replicas replicas summed into: ${effective_bw_for_matrix}."
        fi
        matrix_scale_factor=$(awk "BEGIN {printf \"%.10f\", 1/$num_replicas}")
    fi

    # Step 1: Compute Matrix around TSS
    log "Step 1: Preparing to run computeMatrix for ${sample_group_prefix}."
    local computeMatrix_log="${log_file_dir}/computeMatrix.log"
    local run_compute_matrix=true

    if [[ -f "${matrix_gz}" && -s "${matrix_gz}" ]]; then
        if [[ $num_replicas -eq 1 ]]; then
            log "Found existing matrix file: ${matrix_gz} (for single replica). Will use it."
            run_compute_matrix=false
            # Verify companion files for single replica case
            if [[ ! -f "${matrix_tab}" || ! -s "${matrix_tab}" ]]; then
                log "WARNING: Companion matrix file ${matrix_tab} is missing or empty for existing ${matrix_gz}."
            fi
            if [[ ! -f "${matrix_bed}" || ! -s "${matrix_bed}" ]]; then
                log "WARNING: Companion matrix file ${matrix_bed} is missing or empty for existing ${matrix_gz}."
            fi
        else # num_replicas > 1 (averaging was or should have been involved)
            log "Found existing matrix file: ${matrix_gz}. Since multiple replicas were involved for averaging,"
            log "this matrix will be REGENERATED to ensure it's based on the current single-sample averaging logic."
            rm -f "${matrix_gz}" "${matrix_tab}" "${matrix_bed}" && log "Old matrix files removed." || log "Warning: Failed to remove one or more old matrix files. computeMatrix will attempt to overwrite."
            run_compute_matrix=true
        fi
    else
        log "No existing valid matrix file found (${matrix_gz}), or it is empty. computeMatrix will run."
        run_compute_matrix=true
    fi

    if [[ "$run_compute_matrix" == true ]]; then
        log "Running computeMatrix for ${sample_group_prefix}..."
        log "Effective BigWig for computeMatrix: ${effective_bw_for_matrix}"
        log "Scale factor for computeMatrix: ${matrix_scale_factor}"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b ${COMPUTE_MATRIX_UPSTREAM} \
            -a ${COMPUTE_MATRIX_DOWNSTREAM} \
            -R "${GENE_ANNOTATION_GTF}" \
            -S "${effective_bw_for_matrix}" \
            --scale "${matrix_scale_factor}" \
            --samplesLabel "${sample_group_prefix}" \
            -o "${matrix_gz}" \
            --outFileNameMatrix "${matrix_tab}" \
            --outFileSortedRegions "${matrix_bed}" \
            --skipZeros \
            --missingDataAsZero \
            --binSize ${COMPUTE_MATRIX_BINSIZE} \
            --numberOfProcessors ${COMPUTE_MATRIX_THREADS} \
            --verbose \
            >& "${computeMatrix_log}" \
            || die "computeMatrix failed for ${sample_group_prefix}. Check log: ${computeMatrix_log}"

        [[ -s "$matrix_gz" ]] || die "computeMatrix failed to produce output file: $matrix_gz (or it is empty). Check log: ${computeMatrix_log}"
        log "computeMatrix complete for ${sample_group_prefix}."
    else
        log "Skipped running computeMatrix as a compatible existing matrix was found."
    fi

    # Step 2: Plot Heatmap from the (potentially averaged) Matrix
    log "Step 2: Running plotHeatmap for ${sample_group_prefix}..."
    local plotHeatmap_log="${log_file_dir}/plotHeatmap.log"
    # local PLOT_HEATMAP_TITLE="Signal Heatmap around TSS (${sample_group_prefix} - Averaged Replicas)" # Define title
    local PLOT_HEATMAP_TITLE="TSS ${sample_group_prefix}" # Define 

    # Plot PNG
    # Rely on label embedded in matrix by computeMatrix instead of using --samplesLabel
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
        --plotTitle "${PLOT_HEATMAP_TITLE}" \
        --verbose \
        >& "${plotHeatmap_log}" || die "plotHeatmap (PNG) failed for ${sample_group_prefix}. Check log: ${plotHeatmap_log}"

    # Verify plotHeatmap output
    [[ -s "$heatmap_plot_png" ]] || die "plotHeatmap failed to produce PNG output file: $heatmap_plot_png. Check log: ${plotHeatmap_log}"
    log "plotHeatmap complete for ${sample_group_prefix}."

    # --- Final Summary ---
    log "--- Summary for ${sample_group_prefix} (Merged Replicas) ---"
    log "Input BigWig files used:"
    for bwf in "${input_bw_files[@]}"; do
        log "  - ${bwf}"
    done
    log "Output Directory: ${output_dir}"
    log "Matrix File: ${matrix_gz}"
    log "Heatmap Plot (PNG): ${heatmap_plot_png}"
    log "Temporary files were in: ${SAMPLE_GROUP_TMP_DIR}" # Will be cleaned up by trap
    log "--- Processing Complete for ${sample_group_prefix} ---"

} # End of main function

# --- Execute Main ---
# Make SAMPLE_GROUP_TMP_DIR global and initialize for the trap
SAMPLE_GROUP_TMP_DIR=""
# Pass all script arguments to the main function
main "$@"

# Cleanup is handled by the trap on EXIT/ERR/INT/TERM
log "Script finished successfully."
exit 0