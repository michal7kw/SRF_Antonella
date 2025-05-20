#!/bin/bash
# shellcheck disable=SC2086,SC2046,SC2207

# ==============================================================================
# Script: Generate Comparison Metaprofiles around TSS for GFP and YAF Groups
#
# Description:
#   Calculates and plots the average signal profiles (metaprofiles) of
#   CUT&Tag signal for GFP and YAF sample groups, centered around gene
#   Transcription Start Sites (TSSs). The script processes each group
#   individually to generate a signal matrix (summing replicas if present)
#   and then plots both profiles on a single comparative graph.
#   Uses deepTools computeMatrix and plotProfile.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories for comparison.
#   2. Define input paths (BigWig directory, Gene Annotation).
#   3. For each sample group (GFP, YAF):
#      a. Identify replica BigWig files.
#      b. If multiple replicas, sum them into a single temporary BigWig file.
#      c. Calculate a signal matrix around TSS using computeMatrix with the
#         (potentially summed) BigWig, applying scaling if replicas were summed.
#   4. Generate a single metaprofile plot comparing GFP and YAF using plotProfile
#      with both generated matrices.
#   5. Clean up temporary files.
#
# Usage:
#   bash 13_metaprofile_tss_comparison.sh
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

# Output directory (relative to SCRIPT_DIR) - Adjusted for comparison analysis
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/13_metaprofile_tss/comparison"
LOG_DIR="${SCRIPT_DIR}/logs/13_metaprofile_tss_comparison"
INTERMEDIATES_DIR="${BASE_OUTPUT_DIR}/intermediates" # For persistent intermediate files
TMP_BASE_DIR="${SCRIPT_DIR}/tmp_metaprofile_tss_comparison" # Base for temporary files for this run

# Sample group prefixes to be processed and compared
SAMPLE_GROUP_PREFIXES=(GFP YAF)
# Assign colors for each group for the combined plot.
SAMPLE_GROUP_COLORS=(blue red) # GFP will be blue, YAF will be red. Add more if more groups.

# deepTools Parameters
COMPUTE_MATRIX_UPSTREAM=5000      # Base pairs upstream of TSS
COMPUTE_MATRIX_DOWNSTREAM=5000    # Base pairs downstream of TSS
COMPUTE_MATRIX_BINSIZE=50         # Bin size for averaging signal
COMPUTE_MATRIX_THREADS=40          # Number of threads for computeMatrix
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

    local matrix_files_to_plot=()
    local sample_labels_for_plot=()
    local colors_for_plot_final=() # Renamed to avoid conflict

    # --- Process each sample group (GFP, YAF) to generate its matrix ---
    for i in ${!SAMPLE_GROUP_PREFIXES[@]}; do
        local sample_group_prefix=${SAMPLE_GROUP_PREFIXES[$i]}
        local group_color=${SAMPLE_GROUP_COLORS[$i]:-"gray"} # Default color if not enough

        log "--- Starting matrix generation for sample group: ${sample_group_prefix} ---"

        # Define paths specific to this group's processing
        # Intermediates will be stored persistently
        local group_intermediate_dir="${INTERMEDIATES_DIR}/${sample_group_prefix}"
        local group_log_file_dir="${LOG_DIR}/${sample_group_prefix}"
        # Temporary files for summing steps if not already done, can still use COMPARISON_TMP_DIR
        local group_bw_summing_temp_storage_for_steps="${COMPARISON_TMP_DIR}/${sample_group_prefix}_bw_summing_steps"


        mkdir -p "${group_intermediate_dir}" || die "Failed to create intermediate dir for ${sample_group_prefix}: ${group_intermediate_dir}"
        mkdir -p "${group_log_file_dir}" || die "Failed to create log dir for ${sample_group_prefix}: ${group_log_file_dir}"
        mkdir -p "${group_bw_summing_temp_storage_for_steps}" || die "Failed to create BW summing step temp for ${sample_group_prefix}: ${group_bw_summing_temp_storage_for_steps}"

        # --- Find BigWig files ---
        local input_bw_files_group=() # Use a group-specific array name
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
        local matrix_scale_factor=1.0
        # Path for the final (potentially) summed BigWig for this group - stored in intermediates
        local final_summed_bw_path_group="${group_intermediate_dir}/${sample_group_prefix}_averaged_sum.bw"

        if [[ $num_replicas -eq 0 ]]; then
            die "Logic error: No BigWig files found for ${sample_group_prefix} before processing."
        elif [[ $num_replicas -eq 1 ]]; then
            log "Single replica found for ${sample_group_prefix}. Using it directly: ${input_bw_files_group[0]}"
            effective_bw_for_matrix="${input_bw_files_group[0]}"
        else # num_replicas > 1
            log "Multiple replicas ($num_replicas) found for ${sample_group_prefix}."
            if [[ -f "${final_summed_bw_path_group}" && -s "${final_summed_bw_path_group}" ]]; then
                log "Found existing summed BigWig for ${sample_group_prefix}: ${final_summed_bw_path_group}. Using it."
                effective_bw_for_matrix="${final_summed_bw_path_group}"
            else
                log "No existing summed BigWig for ${sample_group_prefix} found or it is empty. Summing replicas using bigwigCompare..."
                local bigwigCompare_sum_log_group="${group_log_file_dir}/bigwigCompare_sum.log"
                >"${bigwigCompare_sum_log_group}" # Ensure log is fresh

                local sum_in_progress_bw="${input_bw_files_group[0]}"
                # Intermediates for summing steps can use the run-specific temp dir
                mkdir -p "${group_bw_summing_temp_storage_for_steps}"

                for (( k=1; k<num_replicas; k++ )); do
                    local file_to_add="${input_bw_files_group[$k]}"
                    local current_step_output_bw
                    if [[ $k -eq $((num_replicas - 1)) ]]; then # Last iteration, save to final persistent path for this group
                        current_step_output_bw="${final_summed_bw_path_group}"
                    else
                        # Step intermediates are temporary for this summing process
                        current_step_output_bw="${group_bw_summing_temp_storage_for_steps}/${sample_group_prefix}_sum_intermediate_step_${k}.bw"
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
                        >> "${bigwigCompare_sum_log_group}" 2>&1 \
                        || die "bigwigCompare add operation failed for ${sample_group_prefix}. Check log: ${bigwigCompare_sum_log_group}"
                    [[ -s "${current_step_output_bw}" ]] || die "bigwigCompare add operation did not produce output file for ${sample_group_prefix}: ${current_step_output_bw}"
                    
                    sum_in_progress_bw="${current_step_output_bw}"
                done
                effective_bw_for_matrix="${final_summed_bw_path_group}"
                log "All $num_replicas replicas for ${sample_group_prefix} summed into: ${effective_bw_for_matrix}."
            fi
            matrix_scale_factor=$(awk "BEGIN {printf \"%.10f\", 1/$num_replicas}")
        fi

        # --- Compute Matrix for this group ---
        log "Step 1 (Group ${sample_group_prefix}): Preparing to run computeMatrix."
        # Matrix files for this group will be stored in its persistent intermediate directory
        local matrix_gz_group="${group_intermediate_dir}/${sample_group_prefix}_tss_matrix.mat.gz"
        local matrix_tab_group="${group_intermediate_dir}/${sample_group_prefix}_tss_matrix.tab"
        local matrix_bed_group="${group_intermediate_dir}/${sample_group_prefix}_tss_matrix.bed"
        local computeMatrix_log_group="${group_log_file_dir}/computeMatrix.log"
        
        # Decide whether to run computeMatrix or use existing persistent matrix
        local run_compute_matrix_group=true
        if [[ -f "${matrix_gz_group}" && -s "${matrix_gz_group}" ]]; then
            log "Found existing persistent matrix file for ${sample_group_prefix}: ${matrix_gz_group}. Will use it."
            # Also check for companion files, though their absence might not be critical for plotProfile if matrix_gz is fine
            [[ -f "${matrix_tab_group}" ]] || log "WARNING: Companion matrix file ${matrix_tab_group} missing for existing ${matrix_gz_group}."
            [[ -f "${matrix_bed_group}" ]] || log "WARNING: Companion matrix file ${matrix_bed_group} missing for existing ${matrix_gz_group}."
            run_compute_matrix_group=false
        else
            log "No existing valid persistent matrix file found for ${sample_group_prefix} (${matrix_gz_group}), or it is empty. computeMatrix will run."
        fi

        if [[ "$run_compute_matrix_group" == true ]]; then
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
                -o "${matrix_gz_group}" \
                --outFileNameMatrix "${matrix_tab_group}" \
                --outFileSortedRegions "${matrix_bed_group}" \
                --skipZeros \
                --missingDataAsZero \
                --binSize ${COMPUTE_MATRIX_BINSIZE} \
                --numberOfProcessors ${COMPUTE_MATRIX_THREADS} \
                --verbose \
                >& "${computeMatrix_log_group}" \
                || die "computeMatrix failed for ${sample_group_prefix}. Check log: ${computeMatrix_log_group}"

            [[ -s "$matrix_gz_group" ]] || die "computeMatrix failed to produce output file for ${sample_group_prefix}: $matrix_gz_group (or it is empty). Check log: ${computeMatrix_log_group}"
            log "computeMatrix complete for ${sample_group_prefix}."
        fi

        matrix_files_to_plot+=("${matrix_gz_group}")
        sample_labels_for_plot+=("${sample_group_prefix}")
        colors_for_plot_final+=("${group_color}")

        # Optionally clean up group_bw_summing_tmp_dir here if it's large and not needed further
        # For now, the main `cleanup` function will handle the parent `COMPARISON_TMP_DIR`
        log "--- Finished matrix generation for sample group: ${sample_group_prefix} ---"
    done

    # --- Plot Combined Profile ---
    log "Step 2: Running plotProfile for combined comparison of groups: ${sample_labels_for_plot[*]}..."
    # Construct a filename base like "GFP_YAF"
    local filename_base
    filename_base=$(printf "%s_" "${SAMPLE_GROUP_PREFIXES[@]}") # Adds trailing underscore
    filename_base="${filename_base%_}" # Removes trailing underscore
    local final_comparison_plot_png="${BASE_OUTPUT_DIR}/${filename_base}_tss_comparison_profile.png"
    local plotProfile_comparison_log="${LOG_DIR}/plotProfile_comparison.log"
    local plot_title_comparison="TSS Profiles: ${sample_labels_for_plot[*]}"

    # Check if final plot outputs already exist
    if [[ -f "${final_comparison_plot_png}" && -s "${final_comparison_plot_png}" ]]; then
        log "Final comparison profile PNG already exists: ${final_comparison_plot_png}"
        log "Skipping combined plotProfile generation. If you want to regenerate, please delete the existing file."
    else
        log "Plotting combined profile to: ${final_comparison_plot_png}"
        log "Using matrices: ${matrix_files_to_plot[*]}"
        # Labels will be taken from the matrix files directly.
        log "Using colors: ${colors_for_plot_final[*]}"

        # Construct the -m arguments: -m file1 -m file2 ...
        local plot_profile_m_args=()
        for matrix_f in "${matrix_files_to_plot[@]}"; do
            plot_profile_m_args+=("-m" "${matrix_f}")
        done

        plotProfile \
            "${plot_profile_m_args[@]}" \
            --plotType lines \
            -out "${final_comparison_plot_png}" \
            --colors ${colors_for_plot_final[@]} \
            --plotHeight ${PLOT_PROFILE_HEIGHT} \
            --plotWidth ${PLOT_PROFILE_WIDTH} \
            --yAxisLabel "Average Signal" \
            --plotTitle "${plot_title_comparison}" \
            --legendLocation best \
            --verbose \
            >& "${plotProfile_comparison_log}" || die "plotProfile (PNG) for comparison failed. Check log: ${plotProfile_comparison_log}"

        [[ -s "$final_comparison_plot_png" ]] || die "plotProfile failed to produce comparison PNG output file: $final_comparison_plot_png. Check log: ${plotProfile_comparison_log}"
        log "plotProfile for comparison complete. Output: ${final_comparison_plot_png}"
    fi

    # --- Final Summary ---
    log "--- Summary for Comparison Metaprofile (${SAMPLE_GROUP_PREFIXES[*]}) ---"
    log "Input Gene Annotation: ${GENE_ANNOTATION_GTF}"
    log "Input BigWig Directory: ${BIGWIG_DIR}"
    log "Matrices used for plotting (located in ${INTERMEDIATES_DIR}):"
    for matrix_file in "${matrix_files_to_plot[@]}"; do
        log "  - ${matrix_file}"
    done
    log "Output Directory for Plot: ${BASE_OUTPUT_DIR}"
    log "Intermediate files (summed BigWigs, matrices) are in: ${INTERMEDIATES_DIR}"
    log "Comparison Profile Plot (PNG): ${final_comparison_plot_png}"
    log "Log files are in: ${LOG_DIR}"
    log "Run-specific temporary files were in: ${COMPARISON_TMP_DIR} (will be cleaned up on exit)"
    log "--- Comparison Processing Complete ---"

} # End of main function

# --- Execute Main ---
main "$@"

log "Comparison script finished successfully."
exit 0