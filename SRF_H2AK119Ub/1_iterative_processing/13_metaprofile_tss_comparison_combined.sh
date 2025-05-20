#!/bin/bash
# shellcheck disable=SC2086,SC2046,SC2207

# ==============================================================================
# Script: Generate Comparison Metaprofiles around TSS for GFP and YAF Groups
#
# Description:
#   Calculates and plots the average signal profiles (metaprofiles) of
#   CUT&Tag signal for GFP and YAF sample groups, centered around gene
#   Transcription Start Sites (TSSs). This script processes both groups
#   in a single computeMatrix call to ensure proper comparison.
#   Uses deepTools computeMatrix and plotProfile.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Define input paths (BigWig directory, Gene Annotation).
#   3. For each sample group (GFP, YAF):
#      a. Identify replica BigWig files.
#      b. If multiple replicas, sum them into a single BigWig file.
#   4. Run computeMatrix with both summed BigWigs as input.
#   5. Generate a metaprofile plot comparing GFP and YAF.
#   6. Clean up temporary files.
#
# Usage:
#   bash 13_metaprofile_tss_comparison_combined.sh
#   (No arguments needed, processes GFP and YAF by default)
#
# Requirements:
#   - Conda environment with deepTools installed and activated.
#   - Input BigWig files, named like {sample_group_prefix}_{replica_id}.bw
#     (e.g., GFP_1.bw, YAF_1.bw).
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
GENE_ANNOTATION_GTF="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"

# Directory containing the input BigWig files.
# Assumes files are named {sample_group_prefix}_{replica_id}.bw
BIGWIG_DIR="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/10_bigwig"
# === END USER INPUT ===

# Output directory (relative to SCRIPT_DIR)
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/13_metaprofile_tss/comparison_combined"
LOG_DIR="${SCRIPT_DIR}/logs/13_metaprofile_tss_comparison_combined"
INTERMEDIATES_DIR="${BASE_OUTPUT_DIR}/intermediates" # For persistent intermediate files
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_metaprofile_tss_comparison_combined" # Base for temporary files

# Sample group prefixes to be processed and compared
SAMPLE_GROUP_PREFIXES=(GFP YAF)
# Assign colors for each group for the combined plot.
SAMPLE_GROUP_COLORS=(blue red) # GFP will be blue, YAF will be red.

# deepTools Parameters
COMPUTE_MATRIX_UPSTREAM=5000      # Base pairs upstream of TSS
COMPUTE_MATRIX_DOWNSTREAM=5000    # Base pairs downstream of TSS
COMPUTE_MATRIX_BINSIZE=50         # Bin size for averaging signal
COMPUTE_MATRIX_THREADS=40         # Number of threads for computeMatrix
PLOT_PROFILE_HEIGHT=15            # Plot height in cm
PLOT_PROFILE_WIDTH=10             # Plot width in cm (adjusted for two profiles)

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

# This global variable will hold the path to the main temporary directory for this script run.
COMPARISON_TMP_DIR=""

cleanup() {
    # Removes the main temporary directory if it exists
    log "Executing cleanup for comparison script..."
    if [[ -n "${COMPARISON_TMP_DIR:-}" && -d "${COMPARISON_TMP_DIR}" ]]; then
        log "Removing main temporary directory: ${COMPARISON_TMP_DIR}"
        rm -rf "${COMPARISON_TMP_DIR}"
        log "Main temporary directory removed."
    else
         log "Skipping cleanup: Main temporary directory path invalid ('${COMPARISON_TMP_DIR:-}') or does not exist."
    fi
}

# --- Main Script Logic ---

main() {
    # --- Setup Directories for Comparison ---
    COMPARISON_TMP_DIR="${TMP_BASE_DIR}/run_$(date +%Y%m%d%H%M%S%N)" # Assign to global

    mkdir -p "${BASE_OUTPUT_DIR}" || die "Failed to create base output directory: ${BASE_OUTPUT_DIR}"
    mkdir -p "${LOG_DIR}" || die "Failed to create log directory: ${LOG_DIR}"
    mkdir -p "${INTERMEDIATES_DIR}" || die "Failed to create intermediates directory: ${INTERMEDIATES_DIR}"
    mkdir -p "${COMPARISON_TMP_DIR}" || die "Failed to create main temporary directory: ${COMPARISON_TMP_DIR}"
    export TMPDIR="${COMPARISON_TMP_DIR}" # Set for tools that use TMPDIR (e.g. sort)

    # --- Trap for Cleanup ---
    trap cleanup EXIT ERR INT TERM

    log "Starting comparison metaprofile generation for groups: ${SAMPLE_GROUP_PREFIXES[*]}"
    # --- Check Dependencies ---
    check_tool computeMatrix
    check_tool plotProfile
    check_tool bigwigCompare
    check_tool mkdir
    check_tool date
    check_tool rm
    check_tool printf
    check_tool find
    check_tool sort
    check_tool awk

    # --- Verify Input Files/Placeholders ---
    [[ -f "$GENE_ANNOTATION_GTF" ]] || die "Gene annotation file not found: $GENE_ANNOTATION_GTF"
    [[ -d "$BIGWIG_DIR" ]] || die "BigWig directory not found: $BIGWIG_DIR"

    # Arrays to store the final BigWig files for each group
    local final_bigwigs=()
    local sample_labels=()
    local colors=()

    # --- Process each sample group (GFP, YAF) to prepare BigWigs ---
    for i in ${!SAMPLE_GROUP_PREFIXES[@]}; do
        local sample_group_prefix=${SAMPLE_GROUP_PREFIXES[$i]}
        local group_color=${SAMPLE_GROUP_COLORS[$i]:-"gray"} # Default color if not enough

        log "--- Starting BigWig preparation for sample group: ${sample_group_prefix} ---"

        # Define paths specific to this group's processing
        local group_intermediate_dir="${INTERMEDIATES_DIR}/${sample_group_prefix}"
        local group_log_file_dir="${LOG_DIR}/${sample_group_prefix}"
        local group_bw_summing_temp_storage="${COMPARISON_TMP_DIR}/${sample_group_prefix}_bw_summing_steps"

        mkdir -p "${group_intermediate_dir}" || die "Failed to create intermediate dir for ${sample_group_prefix}: ${group_intermediate_dir}"
        mkdir -p "${group_log_file_dir}" || die "Failed to create log dir for ${sample_group_prefix}: ${group_log_file_dir}"
        mkdir -p "${group_bw_summing_temp_storage}" || die "Failed to create BW summing temp for ${sample_group_prefix}: ${group_bw_summing_temp_storage}"

        # --- Find BigWig files ---
        local input_bw_files_group=()
        local found_files_str
        found_files_str=$(find "${BIGWIG_DIR}" -maxdepth 1 -name "${sample_group_prefix}_*.bw" -print0 | sort -zV | xargs -0)

        if [[ -z "$found_files_str" ]]; then
            die "No BigWig files found for prefix '${sample_group_prefix}' in directory '${BIGWIG_DIR}' with pattern '${sample_group_prefix}_*.bw'"
        fi
        read -r -a input_bw_files_group <<< "$found_files_str"

        log "Found ${#input_bw_files_group[@]} BigWig files for sample group ${sample_group_prefix}:"
        for bwf in "${input_bw_files_group[@]}"; do
            log "  - ${bwf}"
            [[ -f "$bwf" ]] || die "Input BigWig file not found: $bwf (listed for group ${sample_group_prefix})"
        done

        # --- Sum replicas if needed ---
        local num_replicas=${#input_bw_files_group[@]}
        local effective_bw_for_matrix
        # Path for the final (potentially) summed BigWig for this group
        local final_summed_bw_path="${group_intermediate_dir}/${sample_group_prefix}_averaged_sum.bw"

        if [[ $num_replicas -eq 0 ]]; then
            die "Logic error: No BigWig files found for ${sample_group_prefix} before processing."
        elif [[ $num_replicas -eq 1 ]]; then
            log "Single replica found for ${sample_group_prefix}. Using it directly: ${input_bw_files_group[0]}"
            effective_bw_for_matrix="${input_bw_files_group[0]}"
        else # num_replicas > 1
            log "Multiple replicas ($num_replicas) found for ${sample_group_prefix}."
            if [[ -f "${final_summed_bw_path}" && -s "${final_summed_bw_path}" ]]; then
                log "Found existing summed BigWig for ${sample_group_prefix}: ${final_summed_bw_path}. Using it."
                effective_bw_for_matrix="${final_summed_bw_path}"
            else
                log "No existing summed BigWig for ${sample_group_prefix} found or it is empty. Summing replicas using bigwigCompare..."
                local bigwigCompare_sum_log="${group_log_file_dir}/bigwigCompare_sum.log"
                >"${bigwigCompare_sum_log}" # Ensure log is fresh

                local sum_in_progress_bw="${input_bw_files_group[0]}"
                mkdir -p "${group_bw_summing_temp_storage}"

                for (( k=1; k<num_replicas; k++ )); do
                    local file_to_add="${input_bw_files_group[$k]}"
                    local current_step_output_bw
                    if [[ $k -eq $((num_replicas - 1)) ]]; then # Last iteration, save to final path
                        current_step_output_bw="${final_summed_bw_path}"
                    else
                        current_step_output_bw="${group_bw_summing_temp_storage}/${sample_group_prefix}_sum_intermediate_step_${k}.bw"
                    fi
                    
                    log "Summing for ${sample_group_prefix} (step $((k+1)) of $num_replicas): ${sum_in_progress_bw} + ${file_to_add} -> ${current_step_output_bw}"
                    bigwigCompare \
                        --bigwig1 "${sum_in_progress_bw}" \
                        --bigwig2 "${file_to_add}" \
                        --operation add \
                        --binSize 10 \
                        --skipNonCoveredRegions \
                        --numberOfProcessors "${COMPUTE_MATRIX_THREADS}" \
                        -o "${current_step_output_bw}" \
                        >> "${bigwigCompare_sum_log}" 2>&1 \
                        || die "bigwigCompare add operation failed for ${sample_group_prefix}. Check log: ${bigwigCompare_sum_log}"
                    [[ -s "${current_step_output_bw}" ]] || die "bigwigCompare add operation did not produce output file for ${sample_group_prefix}: ${current_step_output_bw}"
                    
                    sum_in_progress_bw="${current_step_output_bw}"
                done
                effective_bw_for_matrix="${final_summed_bw_path}"
                log "All $num_replicas replicas for ${sample_group_prefix} summed into: ${effective_bw_for_matrix}."
            fi
        fi

        # Add this group's BigWig to the list for computeMatrix
        final_bigwigs+=("${effective_bw_for_matrix}")
        sample_labels+=("${sample_group_prefix}")
        colors+=("${group_color}")

        log "--- Finished BigWig preparation for sample group: ${sample_group_prefix} ---"
    done

    # --- Run computeMatrix with all BigWigs ---
    log "Running computeMatrix for all groups: ${sample_labels[*]}..."
    
    # Construct a filename base like "GFP_YAF"
    local filename_base
    filename_base=$(printf "%s_" "${SAMPLE_GROUP_PREFIXES[@]}") # Adds trailing underscore
    filename_base="${filename_base%_}" # Removes trailing underscore
    
    # Output files
    local matrix_gz="${BASE_OUTPUT_DIR}/${filename_base}_tss_matrix.mat.gz"
    local matrix_tab="${BASE_OUTPUT_DIR}/${filename_base}_tss_matrix.tab"
    local matrix_bed="${BASE_OUTPUT_DIR}/${filename_base}_tss_matrix.bed"
    local final_plot_png="${BASE_OUTPUT_DIR}/${filename_base}_tss_comparison_profile.png"
    local computeMatrix_log="${LOG_DIR}/computeMatrix.log"
    local plotProfile_log="${LOG_DIR}/plotProfile.log"
    
    # Check if final outputs already exist
    if [[ -f "${final_plot_png}" && -s "${final_plot_png}" ]]; then
        log "Final plot already exists: ${final_plot_png}"
        log "To regenerate, please delete the existing file."
        exit 0
    fi
    
    # Construct the -S arguments for computeMatrix: -S file1 -S file2 ...
    local compute_matrix_s_args=()
    for bw_file in "${final_bigwigs[@]}"; do
        compute_matrix_s_args+=("-S" "${bw_file}")
    done
    
    # Run computeMatrix
    log "Running computeMatrix with BigWigs: ${final_bigwigs[*]}"
    log "Using sample labels: ${sample_labels[*]}"

    # Build the computeMatrix command parts
    local cmd_args=(reference-point \
        --referencePoint TSS \
        -b "${COMPUTE_MATRIX_UPSTREAM}" \
        -a "${COMPUTE_MATRIX_DOWNSTREAM}" \
        -R "${GENE_ANNOTATION_GTF}")

    # Add score files (-S arguments)
    cmd_args+=(${compute_matrix_s_args[@]})

    # Add sample labels
    cmd_args+=(--samplesLabel)
    for label in "${sample_labels[@]}"; do
        cmd_args+=("${label}")
    done

    # Add output and other parameters
    cmd_args+=(\
        -o "${matrix_gz}" \
        --outFileNameMatrix "${matrix_tab}" \
        --outFileSortedRegions "${matrix_bed}" \
        --skipZeros \
        --missingDataAsZero \
        --binSize "${COMPUTE_MATRIX_BINSIZE}" \
        --numberOfProcessors "${COMPUTE_MATRIX_THREADS}" \
        --verbose)
    
    log "Executing computeMatrix with arguments: ${cmd_args[*]}"
    computeMatrix "${cmd_args[@]}" >& "${computeMatrix_log}" \
        || die "computeMatrix failed. Check log: ${computeMatrix_log}"
    
    [[ -s "$matrix_gz" ]] || die "computeMatrix failed to produce output file: $matrix_gz (or it is empty). Check log: ${computeMatrix_log}"
    log "computeMatrix complete."
    
    # --- Plot Profile ---
    log "Running plotProfile..."
    local plot_title="TSS Profiles: ${sample_labels[*]}"
    
    # Construct the --colors arguments: --colors color1 color2 ...
    local colors_args=()
    for color in "${colors[@]}"; do
        colors_args+=("${color}")
    done
    
    plotProfile \
        -m "${matrix_gz}" \
        --plotType lines \
        -out "${final_plot_png}" \
        --colors ${colors_args[@]} \
        --plotHeight ${PLOT_PROFILE_HEIGHT} \
        --plotWidth ${PLOT_PROFILE_WIDTH} \
        --yAxisLabel "Average Signal" \
        --plotTitle "${plot_title}" \
        --legendLocation best \
        --verbose \
        >& "${plotProfile_log}" || die "plotProfile failed. Check log: ${plotProfile_log}"
    
    [[ -s "$final_plot_png" ]] || die "plotProfile failed to produce output file: $final_plot_png. Check log: ${plotProfile_log}"
    log "plotProfile complete. Output: ${final_plot_png}"
    
    # --- Final Summary ---
    log "--- Summary for Comparison Metaprofile (${SAMPLE_GROUP_PREFIXES[*]}) ---"
    log "Input Gene Annotation: ${GENE_ANNOTATION_GTF}"
    log "Input BigWig Directory: ${BIGWIG_DIR}"
    log "BigWigs used for plotting:"
    for i in "${!final_bigwigs[@]}"; do
        log "  - ${sample_labels[$i]}: ${final_bigwigs[$i]} (${colors[$i]})"
    done
    log "Output Directory: ${BASE_OUTPUT_DIR}"
    log "Matrix File: ${matrix_gz}"
    log "Profile Plot: ${final_plot_png}"
    log "Log files are in: ${LOG_DIR}"
    log "--- Comparison Processing Complete ---"

} # End of main function

# --- Execute Main ---
main "$@"

log "Comparison script finished successfully."
exit 0