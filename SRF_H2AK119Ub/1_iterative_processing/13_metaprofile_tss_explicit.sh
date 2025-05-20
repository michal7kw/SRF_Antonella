#!/bin/bash
# shellcheck disable=SC2086,SC2046,SC2207

# ==============================================================================
# Script: Generate Comparison Metaprofiles around TSS for GFP and YAF Groups (Averaging Replicates)
#
# Description:
#   Calculates and plots the average signal profiles (metaprofiles) of
#   CUT&Tag signal for GFP and YAF sample groups, centered around gene
#   Transcription Start Sites (TSSs). If multiple replicates exist for a
#   sample group, their signals are summed and averaged.
#   Uses deepTools computeMatrix, plotProfile, and bigwigCompare.
#
# Usage:
#   bash 13_metaprofile_tss_explicit.sh
# ==============================================================================

# --- Configuration ---
set -euo pipefail

# Script directory for relative paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- Parameters & Paths ---
GENE_ANNOTATION_GTF="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"
BIGWIG_DIR="/mnt/d/Github/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/10_bigwig"

# Output directory
OUTPUT_DIR="${SCRIPT_DIR}/analysis/13_metaprofile_tss/explicit"
LOG_DIR="${SCRIPT_DIR}/logs/13_metaprofile_tss_explicit"
TMP_WORK_DIR="${SCRIPT_DIR}/tmp_13_metaprofile_tss_explicit_work" # Temporary working directory

# deepTools Parameters
UPSTREAM=5000
DOWNSTREAM=5000
BINSIZE=50
THREADS=40 # Used for computeMatrix and bigwigCompare
PLOT_HEIGHT=15
PLOT_WIDTH=10

# --- Helper Functions ---
log() {
    printf "[%(%Y-%m-%d %H:%M:%S)T] %s\n" -1 "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

check_tool() {
    command -v "$1" >/dev/null 2>&1 || die "Required tool '$1' not found in PATH. Ensure it is installed and accessible."
}

cleanup() {
    log "Executing cleanup..."
    if [[ -d "${TMP_WORK_DIR}" ]]; then
        log "Removing temporary directory: ${TMP_WORK_DIR}"
        rm -rf "${TMP_WORK_DIR}"
        log "Temporary directory removed."
    else
        log "Skipping cleanup: Temporary directory '${TMP_WORK_DIR}' does not exist."
    fi
}

# --- Helper function to process a sample group's BigWigs ---
# Arguments:
#   $1: Sample group prefix (e.g., "GFP")
#   $2: String containing list of BigWig files (newline-separated)
#   $3: Temporary directory path for summed files
#   $4: Log directory for bigwigCompare logs
#   $5: Number of threads for bigwigCompare
# Outputs:
#   Echoes two values: effective_bw_path scale_factor
process_sample_group() {
    local group_prefix="$1"
    local bigwig_list_str="$2"
    local tmp_dir_summed_files="$3" # Specific tmp dir for this group's summed file
    local log_dir_bwcompare="$4"    # Log dir for this group's bigwigCompare log
    local threads_bwcompare="$5"

    local -a bw_array
    # Read newline-separated string into array, handle empty string if no files found.
    if [[ -n "$bigwig_list_str" ]]; then
        # shellcheck disable=SC2207
        bw_array=($(echo "$bigwig_list_str" | tr '\n' ' '))
    else
        bw_array=()
    fi
    local num_replicas=${#bw_array[@]}
    local effective_bw_path
    local scale_factor

    log "Processing sample group: ${group_prefix} with ${num_replicas} replica(s) found initially."

    if [[ $num_replicas -eq 0 ]]; then
        die "No BigWig files found for group '${group_prefix}' with pattern '${group_prefix}_*.bw' in ${BIGWIG_DIR}"
    fi

    for bwf_idx in "${!bw_array[@]}"; do
        local bwf="${bw_array[$bwf_idx]}"
        log "  - Verifying replica for ${group_prefix}: ${bwf}"
        [[ -f "$bwf" ]] || die "Input BigWig file not found: $bwf (listed for group ${group_prefix})"
    done

    if [[ $num_replicas -eq 1 ]]; then
        log "Single replica found for ${group_prefix}. Using it directly."
        effective_bw_path="${bw_array[0]}"
        scale_factor="1.0"
    else # num_replicas > 1
        log "Multiple replicas ($num_replicas) found for ${group_prefix}. Summing and averaging..."
        local final_summed_bw_path="${tmp_dir_summed_files}/${group_prefix}_averaged_sum.bw"
        local bigwigCompare_sum_log="${log_dir_bwcompare}/bigwigCompare_sum_${group_prefix}.log"
        mkdir -p "$(dirname "${final_summed_bw_path}")" # Ensure tmp_dir_summed_files exists
        >"${bigwigCompare_sum_log}" # Ensure log is fresh

        # Check if already summed (e.g. for reruns, though less likely in this script's context)
        if [[ -f "${final_summed_bw_path}" && -s "${final_summed_bw_path}" ]]; then
            log "Found existing summed BigWig for ${group_prefix}: ${final_summed_bw_path}. Using it."
            effective_bw_path="${final_summed_bw_path}"
        else
            log "Summing ${num_replicas} replicas for ${group_prefix} using bigwigCompare..."
            local sum_in_progress_bw="${bw_array[0]}"
            # Create a subdirectory within tmp_dir_summed_files for intermediate sums to avoid clutter
            local temp_intermediate_sum_storage="${tmp_dir_summed_files}/sum_intermediates_${group_prefix}"
            mkdir -p "${temp_intermediate_sum_storage}"

            for (( i=1; i<num_replicas; i++ )); do
                local file_to_add="${bw_array[$i]}"
                local current_step_output_bw
                if [[ $i -eq $((num_replicas - 1)) ]]; then # Last iteration, save to final path
                    current_step_output_bw="${final_summed_bw_path}"
                else
                    current_step_output_bw="${temp_intermediate_sum_storage}/${group_prefix}_sum_intermediate_step_${i}.bw"
                fi
                
                log "Summing (step $((i+1)) of ${num_replicas} for ${group_prefix}): ${sum_in_progress_bw} + ${file_to_add} -> ${current_step_output_bw}"
                bigwigCompare \
                    --bigwig1 "${sum_in_progress_bw}" \
                    --bigwig2 "${file_to_add}" \
                    --operation add \
                    --binSize 10 \
                    --skipNonCoveredRegions \
                    --numberOfProcessors "${threads_bwcompare}" \
                    -o "${current_step_output_bw}" \
                    >> "${bigwigCompare_sum_log}" 2>&1 \
                    || die "bigwigCompare add operation failed for ${group_prefix}. Check log: ${bigwigCompare_sum_log}"
                [[ -s "${current_step_output_bw}" ]] || die "bigwigCompare add operation did not produce output file for ${group_prefix}: ${current_step_output_bw}"
                
                sum_in_progress_bw="${current_step_output_bw}"
            done
            effective_bw_path="${final_summed_bw_path}"
            log "All $num_replicas replicas for ${group_prefix} summed into: ${effective_bw_path}."
        fi
        scale_factor=$(awk "BEGIN {printf \"%.10f\", 1/$num_replicas}")
    fi
    log "Effective BigWig for ${group_prefix}: ${effective_bw_path}"
    log "Scale factor for ${group_prefix}: ${scale_factor}"
    echo "${effective_bw_path} ${scale_factor}"
}


# --- Main Script ---
main() {
    trap cleanup EXIT ERR INT TERM

    # Check for required tools
    check_tool computeMatrix
    check_tool plotProfile
    check_tool bigwigCompare
    check_tool awk
    check_tool find
    check_tool sort
    check_tool echo
    check_tool tr
    check_tool head # Not used anymore but good to have if needed
    check_tool rm
    check_tool mkdir

    # Create directories
    mkdir -p "${OUTPUT_DIR}" || die "Failed to create output directory"
    mkdir -p "${LOG_DIR}" || die "Failed to create log directory"
    mkdir -p "${TMP_WORK_DIR}" || die "Failed to create temporary working directory"
    
    # Find BigWig files for GFP and YAF
    GFP_BIGWIGS_STR=$(find "${BIGWIG_DIR}" -maxdepth 1 -name "GFP_*.bw" -print0 | sort -zV | xargs -0 printf "%s\n" || true)
    YAF_BIGWIGS_STR=$(find "${BIGWIG_DIR}" -maxdepth 1 -name "YAF_*.bw" -print0 | sort -zV | xargs -0 printf "%s\n" || true)
    
    # Process GFP samples
    local GFP_EFFECTIVE_BW GFP_SCALE_FACTOR
    read -r GFP_EFFECTIVE_BW GFP_SCALE_FACTOR <<< "$(process_sample_group "GFP" "$GFP_BIGWIGS_STR" "${TMP_WORK_DIR}/GFP" "$LOG_DIR" "$THREADS")"
    [[ -n "$GFP_EFFECTIVE_BW" && -n "$GFP_SCALE_FACTOR" ]] || die "Failed to process GFP samples."

    # Process YAF samples
    local YAF_EFFECTIVE_BW YAF_SCALE_FACTOR
    read -r YAF_EFFECTIVE_BW YAF_SCALE_FACTOR <<< "$(process_sample_group "YAF" "$YAF_BIGWIGS_STR" "${TMP_WORK_DIR}/YAF" "$LOG_DIR" "$THREADS")"
    [[ -n "$YAF_EFFECTIVE_BW" && -n "$YAF_SCALE_FACTOR" ]] || die "Failed to process YAF samples."

    # Output files
    MATRIX_FILE="${OUTPUT_DIR}/GFP_YAF_tss_matrix.mat.gz"
    PLOT_FILE="${OUTPUT_DIR}/GFP_YAF_tss_comparison.png"
    
    # Run computeMatrix with potentially averaged files and scale factors
    log "Running computeMatrix..."
    computeMatrix reference-point \
        --referencePoint TSS \
        -b $UPSTREAM \
        -a $DOWNSTREAM \
        -R "$GENE_ANNOTATION_GTF" \
        -S "$GFP_EFFECTIVE_BW" "$YAF_EFFECTIVE_BW" \
        --scale "$GFP_SCALE_FACTOR" "$YAF_SCALE_FACTOR" \
        --samplesLabel "GFP_avg" "YAF_avg" \
        -o "$MATRIX_FILE" \
        --skipZeros \
        --missingDataAsZero \
        --binSize $BINSIZE \
        --numberOfProcessors $THREADS \
        --verbose \
        >& "${LOG_DIR}/computeMatrix.log" || die "computeMatrix failed"
    
    # Run plotProfile
    log "Running plotProfile..."
    plotProfile \
        -m "$MATRIX_FILE" \
        --plotType lines \
        -out "$PLOT_FILE" \
        --colors blue red \
        --plotHeight $PLOT_HEIGHT \
        --plotWidth $PLOT_WIDTH \
        --yAxisLabel "Average Signal" \
        --plotTitle "TSS Profiles: GFP (avg) vs YAF (avg)" \
        --legendLocation best \
        --verbose \
        >& "${LOG_DIR}/plotProfile.log" || die "plotProfile failed"
    
    log "Done! Output files:"
    log "  - Matrix: $MATRIX_FILE"
    log "  - Plot: $PLOT_FILE"
    log "Temporary files were in: ${TMP_WORK_DIR}"
}

main "$@"
exit 0