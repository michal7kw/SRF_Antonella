#!/bin/bash
# shellcheck disable=SC2086,SC2046

# ==============================================================================
# Script: Peak Calling (Strict) for CUT&Tag Data
#
# Description:
#   Performs strict broad peak calling and filtering on CUT&Tag data using MACS2.
#   Processes one sample at a time based on a command-line index.
#
# Pipeline Steps:
#   1. Setup environment, variables, and temporary directories.
#   2. Standardize chromosome names in the input BAM file.
#   3. Call broad peaks using MACS2 with strict parameters.
#   4. Filter raw peaks based on fold enrichment and length.
#   5. Validate filtered BED file format.
#   6. Remove peaks overlapping blacklisted regions.
#   7. Generate a visualization-friendly BED file (normalized scores).
#   8. Calculate and report QC metrics (FRiP, peak count, etc.).
#   9. Generate peak statistics file.
#  10. Clean up temporary files.
#
# Usage:
#   bash 5_peak_calling_strict_v2.sh <sample_index>
# Example:
#   bash 5_peak_calling_strict_v2.sh 0  # Processes the first sample (GFP_1)
#
# Parallel Execution:
#   Use GNU Parallel or xargs for parallel processing across samples:
#   N=6 # Number of parallel jobs
#   seq 0 5 | parallel -j $N bash 5_peak_calling_strict_v2.sh {}
#   # OR
#   # printf "%s\n" {0..5} | xargs -I {} -P $N bash 5_peak_calling_strict_v2.sh {}
#
# Requirements:
#   - Conda environment with MACS2, samtools, bedtools installed and activated.
#   - Input BAM files (deduplicated) and BAI index files.
#   - Genome blacklist file (BED format).
# ==============================================================================

# --- Configuration ---
# Exit on error, treat unset variables as errors, propagate exit status through pipes
set -euo pipefail

# Script directory for relative paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- Parameters & Paths (Modify as needed) ---
# Input data directory structure (relative to SCRIPT_DIR, assuming it's inside 1_iterative_processing)
ALIGNMENT_DIR="${SCRIPT_DIR}/analysis/3_alignment"
# Output directory (relative to SCRIPT_DIR, assuming it's inside 1_iterative_processing)
BASE_OUTPUT_DIR="${SCRIPT_DIR}/analysis/5_peak_calling_strict_v2"
LOG_DIR="${SCRIPT_DIR}/logs/5_peak_calling_strict_v2"
TMP_BASE_DIR="/tmp/srf_h2ak_strict_peak_processing_tmp" # Base for temporary files, moved to Linux /tmp for better I/O

# Blacklist file path (relative to SCRIPT_DIR or absolute)
# Assuming COMMON_DATA is two levels up from 1_iterative_processing
BLACKLIST_BED="${SCRIPT_DIR}/../../COMMON_DATA/hg38-blacklist.v2.bed"

# Sample names array
SAMPLES=(GFP_1 GFP_2 GFP_3 YAF_1 YAF_2 YAF_3)

# MACS2 Parameters (Strict)
MACS2_GENOME="hs"
MACS2_FORMAT="BAMPE"
MACS2_QVALUE="0.01"
MACS2_BROAD_CUTOFF="0.01"
MACS2_MIN_LEN="1000" # Minimum peak length
MACS2_KEEP_DUP="1" # Keep duplicates (appropriate for CUT&Tag)
MACS2_BUFFER_SIZE="1000000"
MACS2_SLOCAL="1000"
MACS2_LLOCAL="10000"

# Filtering Parameters (Strict)
MIN_FOLD_ENRICHMENT=5
MIN_PEAK_LEN=1000
MAX_PEAK_LEN=50000

# Visualization Parameters
VIZ_MAX_FOLD=50 # Cap fold enrichment for scaling visualization score

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
    check_tool macs2
    check_tool samtools
    check_tool bedtools
    check_tool awk
    check_tool sed
    check_tool date
    check_tool mkdir
    check_tool wc
    check_tool mv
    check_tool cat
    check_tool rm
    check_tool printf # Used for logging/output

    # --- Setup Directories & Paths ---
    local output_dir="${BASE_OUTPUT_DIR}"
    local log_dir="${LOG_DIR}/${sample_name}" # Sample-specific log dir (optional)
    local timestamp=$(date +%Y%m%d%H%M%S%N) # Added nanoseconds for higher uniqueness
    # Define SAMPLE_TMP_DIR globally for trap cleanup
    # Ensure TMP_BASE_DIR exists before creating sample temp dir
    mkdir -p "${TMP_BASE_DIR}" || die "Failed to create base temporary directory: ${TMP_BASE_DIR}"
    SAMPLE_TMP_DIR="${TMP_BASE_DIR}/${sample_name}_${timestamp}"

    mkdir -p "${output_dir}" || die "Failed to create output directory: ${output_dir}"
    mkdir -p "${log_dir}" || die "Failed to create log directory: ${log_dir}"
    mkdir -p "${SAMPLE_TMP_DIR}" || die "Failed to create sample temporary directory: ${SAMPLE_TMP_DIR}"
    # Set TMPDIR for MACS2 internal use if needed, though --tempdir is preferred
    export TMPDIR="${SAMPLE_TMP_DIR}"

    # --- Trap for Cleanup ---
    # Ensure SAMPLE_TMP_DIR is cleaned up on exit, error, interrupt, or termination
    # Needs SAMPLE_TMP_DIR to be defined *before* trap is set.
    trap cleanup EXIT ERR INT TERM

    # --- Define File Paths ---
    local input_bam="${ALIGNMENT_DIR}/${sample_name}.dedup.bam"
    local input_bai="${ALIGNMENT_DIR}/${sample_name}.dedup.bam.bai"
    local chr_bam="${SAMPLE_TMP_DIR}/${sample_name}.chr.bam"
    local chr_header="${SAMPLE_TMP_DIR}/${sample_name}_header.sam"
    local raw_header_intermediate="${SAMPLE_TMP_DIR}/${sample_name}_raw_header_intermediate.sam" # For decoupling samtools and sed
    local macs2_tmp="${SAMPLE_TMP_DIR}/macs2_tmp_${RANDOM}" # MACS2 specific temp within sample temp

    local raw_peaks="${output_dir}/${sample_name}_broad_peaks.broadPeak"
    local fixed_chrom_bed="${SAMPLE_TMP_DIR}/fixed_chrom.bed"
    local filtered_peaks="${output_dir}/${sample_name}_broad_peaks_filtered.broadPeak"
    local final_peaks="${output_dir}/${sample_name}_broad_peaks_final.broadPeak"
    local viz_peaks="${output_dir}/${sample_name}_broad_peaks_viz.broadPeak"
    local metrics_csv="${output_dir}/${sample_name}_metrics.csv"
    local stats_txt="${output_dir}/${sample_name}_peak_stats.txt"

    # --- Verify Input Files ---
    [[ -f "$input_bam" ]] || die "Input BAM not found: $input_bam"
    [[ -f "$input_bai" ]] || die "Input BAI not found: $input_bai"
    [[ -f "$BLACKLIST_BED" ]] || die "Blacklist BED not found: $BLACKLIST_BED"

    # --- Configure Python Environment for MACS2 ---
    # These might help with memory issues on some systems for MACS2
    log "Configuring Python environment variables for MACS2..."
    export PYTHONPATH=""
    export MACS2_LARGE_DATA=1
    export PYTHONMALLOC=malloc
    export PYTHONWARNINGS=ignore
    export PYTHONHASHSEED=0
    export PYTHONDONTWRITEBYTECODE=1
    unset PYTHONSTARTUP PYTHONHOME PYTHONEXECUTABLE # Clean potential conflicts

    # --- Pipeline Steps ---

    # Step 1: Standardize Chromosome Names (Add 'chr' prefix)
    log "Step 1: Standardizing chromosome names for ${sample_name}..."
    # Step 1a: Dump raw BAM header to an intermediate file
    log "Step 1a: Dumping raw BAM header for ${sample_name} to ${raw_header_intermediate}..."
    samtools view -H "$input_bam" > "$raw_header_intermediate" \
        || die "Failed to dump raw BAM header for $sample_name to $raw_header_intermediate. samtools exit code: $?"

    # Step 1b: Process the header using sed from the intermediate file
    log "Step 1b: Processing header with sed for ${sample_name} from ${raw_header_intermediate}..."
    sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' "$raw_header_intermediate" > "$chr_header" \
        || die "Failed to process header with sed for $sample_name from $raw_header_intermediate. sed exit code: $?"

    # Clean up the intermediate header file
    rm "$raw_header_intermediate"
    samtools reheader "$chr_header" "$input_bam" > "$chr_bam" || die "Failed to reheader BAM for $sample_name"
    samtools index "$chr_bam" || die "Failed to index chr-prefixed BAM for $sample_name"
    log "Chromosome standardization complete for ${sample_name}."

    # Step 2: MACS2 Peak Calling (Strict Broad)
    log "Step 2: Running MACS2 peak calling for ${sample_name}..."
    mkdir -p "$macs2_tmp" || die "Failed to create MACS2 temp dir: $macs2_tmp"
    # Run MACS2 and capture stderr to a log file
    local macs2_log="${log_dir}/macs2_callpeak.log"
    macs2 callpeak \
        -t "$chr_bam" \
        -f "$MACS2_FORMAT" \
        -g "$MACS2_GENOME" \
        --broad \
        -n "${sample_name}_broad" \
        --outdir "$output_dir" \
        -q "$MACS2_QVALUE" \
        --broad-cutoff "$MACS2_BROAD_CUTOFF" \
        --keep-dup "$MACS2_KEEP_DUP" \
        --min-length "$MACS2_MIN_LEN" \
        --bdg \
        --buffer-size "$MACS2_BUFFER_SIZE" \
        --verbose 3 \
        --tempdir "$macs2_tmp" \
        --nomodel \
        --nolambda \
        --slocal "$MACS2_SLOCAL" \
        --llocal "$MACS2_LLOCAL" \
        >& "$macs2_log" \
        || die "MACS2 peak calling failed for ${sample_name}. Check log: $macs2_log"

    # Verify MACS2 output
    [[ -s "$raw_peaks" ]] || die "MACS2 failed to produce output file: $raw_peaks. Check log: $macs2_log"
    log "MACS2 peak calling complete for ${sample_name}."

    # Step 3: Fix Chromosome Names in Peak File (if MACS2 didn't use 'chr')
    log "Step 3: Fixing chromosome names in raw peak file (${raw_peaks})..."
    awk 'BEGIN{OFS="\t"} {
        if ($1 !~ /^chr/) {
            if ($1 ~ /^[0-9XY]+$/) { $1 = "chr" $1 }
            else if ($1 == "MT") { $1 = "chrM" }
            else if ($1 ~ /^(GL|KI|JH|ML|NC_|NW_|NT_)/) { next } # Skip common contigs/patches/unplaced
            else { print "Warning: Skipping unexpected chromosome format (pre-chr check): " $1 " at line " NR > "/dev/stderr"; next }
        }
        # At this point, $1 should start with "chr" or have been skipped.
        # Add a check for $1 being *only* "chr"
        if ($1 == "chr") {
            print "Warning: Skipping line with standalone '\''chr'\'' chromosome: " $1 " at line " NR > "/dev/stderr";
            next;
        }
        print $0
    }' "$raw_peaks" > "$fixed_chrom_bed" || die "Failed to fix chromosome names in peak file for $sample_name"
    log "Chromosome name fixing complete for ${sample_name}."

    # Step 4: Filter Peaks (Strict Criteria) and Standardize BED6+3 Format
    log "Step 4: Filtering peaks (Fold>=${MIN_FOLD_ENRICHMENT}, Len ${MIN_PEAK_LEN}-${MAX_PEAK_LEN}bp) for ${sample_name}..."
    awk -v sample="$sample_name" \
        -v min_fold="$MIN_FOLD_ENRICHMENT" \
        -v min_len="$MIN_PEAK_LEN" \
        -v max_len="$MAX_PEAK_LEN" \
        'BEGIN{OFS="\t"; peak_idx=0} {
        len = $3 - $2;
        orig_fold = $7;
        orig_pval = $8; # -log10(pvalue)
        orig_qval = $9; # -log10(qvalue)

        # Apply filtering criteria
        if(orig_fold >= min_fold && len >= min_len && len <= max_len) {
            peak_idx++;
            # Create standardized BED name
            name = sprintf("%s_peak_%d", sample, peak_idx);
            gsub(/[^a-zA-Z0-9_.-]/, "_", name); # Sanitize name further (allow dot, hyphen)

            # Standardize score (0-1000) - Use q-value for score
            # q-value is -log10(q), so higher is better. Scale linearly.
            # Cap at e.g. -log10(1e-100) = 100 for scaling purposes if needed, but simple scaling is often fine.
            score = int(orig_qval * 10);
            if(score > 1000) score = 1000;
            if(score < 0) score = 0; # Should not happen with -log10 values

            # Standardize strand
            strand = ($6 == "+" || $6 == "-") ? $6 : ".";

            # Output BED6 + 3 standard columns (signalValue, pValue, qValue)
            # Using original fold, pval, qval for these extra columns
            print $1, $2, $3, name, score, strand, orig_fold, orig_pval, orig_qval;
        }
    }' "$fixed_chrom_bed" > "$filtered_peaks" || die "Failed to filter peaks for $sample_name"
    log "Peak filtering complete for ${sample_name}."

    # Step 5: Validate Filtered BED Format
    log "Step 5: Validating filtered BED format (${filtered_peaks})..."
    local filtered_peaks_tmp="${filtered_peaks}.tmp"
    awk 'BEGIN{OFS="\t"; errors=0} {
        if(NF < 6) {
            print "BED Format Error (Columns): Less than 6 fields at line " NR ": " $0 > "/dev/stderr"; errors++; next
        }
        if($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || int($2) < 0 || int($3) < int($2)) {
            print "BED Format Error (Coords): Invalid coordinates at line " NR ": " $0 > "/dev/stderr"; errors++; next
        }
        if($5 !~ /^[0-9]+$/ || int($5) < 0 || int($5) > 1000) {
            print "BED Format Error (Score): Invalid score (col 5) at line " NR ": " $0 > "/dev/stderr"; errors++; next
        }
        if($6 != "." && $6 != "+" && $6 != "-") {
             print "BED Format Error (Strand): Invalid strand (col 6) at line " NR ": " $0 > "/dev/stderr"; errors++; next
        }
        # Optional: Check extra columns if they exist (e.g., numeric fold, pval, qval)
        if (NF >= 9) {
             if ($7 !~ /^[0-9.-]+([eE][+-]?[0-9]+)?$/) { print "BED Format Warning (Col 7): Non-numeric signalValue at line " NR ": " $7 > "/dev/stderr"; }
             if ($8 !~ /^[0-9.-]+([eE][+-]?[0-9]+)?$/) { print "BED Format Warning (Col 8): Non-numeric pValue at line " NR ": " $8 > "/dev/stderr"; }
             if ($9 !~ /^[0-9.-]+([eE][+-]?[0-9]+)?$/) { print "BED Format Warning (Col 9): Non-numeric qValue at line " NR ": " $9 > "/dev/stderr"; }
        }
        print $0
    } END { if(errors>0) { exit 1 } }' "$filtered_peaks" > "$filtered_peaks_tmp" \
        && mv "$filtered_peaks_tmp" "$filtered_peaks" \
        || die "Filtered BED file validation failed for ${sample_name}. Check stderr for details."
    log "BED validation complete for ${sample_name}."

    # Step 6: Remove Blacklisted Regions
    log "Step 6: Removing blacklisted regions from ${filtered_peaks}..."
    local blacklist_log="${log_dir}/bedtools_intersect_blacklist.log"
    bedtools intersect -v -a "$filtered_peaks" -b "$BLACKLIST_BED" -nonamecheck > "$final_peaks" 2> "$blacklist_log" \
        || die "bedtools intersect failed while removing blacklist regions for ${sample_name}. Check log: $blacklist_log"
    # Check if bedtools actually produced output, even if exit code was 0
    # [[ -s "$final_peaks" ]] || log "Warning: Final peak file is empty after blacklist removal for ${sample_name}." # This might be expected
    log "Blacklist removal complete for ${sample_name}."

    # Step 7: Create Visualization File (BED6+3 with Scaled Score)
    log "Step 7: Creating visualization file (${viz_peaks})..."
    awk -v max_fold_viz="$VIZ_MAX_FOLD" 'BEGIN{OFS="\t"} {
        # Scale score based on fold enrichment for visualization (0-1000)
        # Use original fold enrichment (column 7) for scaling
        fold_val = $7;
        scaled_fold = (fold_val > max_fold_viz) ? max_fold_viz : fold_val;
        # Avoid division by zero if max_fold_viz is 0, though it should be > 0
        if (max_fold_viz <= 0) { viz_score = 0 }
        else { viz_score = int((scaled_fold / max_fold_viz) * 1000); }

        if(viz_score > 1000) viz_score = 1000;
        if(viz_score < 0) viz_score = 0;

        # Output standard BED6 + 3 (signalValue, pValue, qValue from original filtered file)
        # Use the calculated viz_score for the 5th column (score)
        print $1, $2, $3, $4, viz_score, $6, $7, $8, $9;
    }' "$final_peaks" > "$viz_peaks" || die "Failed to create visualization peak file for $sample_name"
    log "Visualization file created for ${sample_name}."

    # Step 8: Calculate QC Metrics
    log "Step 8: Calculating QC metrics for ${sample_name}..."
    local total_reads="NA"
    local reads_in_peaks="NA"
    local frip="NA"
    local peak_count=0
    local mean_length=0
    local mean_fold=0

    # Get total reads from original BAM
    total_reads=$(samtools view -c "$input_bam") || { log "Warning: Failed to get total read count from $input_bam"; total_reads="NA"; }

    # Calculate FRiP using chr-prefixed BAM and final peaks
    if [[ -f "$chr_bam" && -s "$final_peaks" ]]; then
        log "Calculating reads in peaks using bedtools intersect..."
        local frip_log="${log_dir}/bedtools_intersect_frip.log"
        # Count read pairs where *both* ends overlap the same peak by at least 50% (-f 0.5 -F 0.5 -e)
        # Note: -u reports each read pair once if it overlaps *any* peak.
        #       If using single-end reads, adjust intersect options.
        reads_in_peaks=$(bedtools intersect -u -a "$chr_bam" -b "$final_peaks" -f 0.5 -F 0.5 -e 2> "$frip_log" | samtools view -c -)

        if [[ $? -ne 0 ]]; then
             log "Warning: bedtools intersect failed during FRiP calculation for ${sample_name}. Check log: $frip_log"
             reads_in_peaks="NA"
             frip="NA"
        elif [[ "$total_reads" != "NA" && "$total_reads" -gt 0 ]]; then
            frip=$(awk -v r="$reads_in_peaks" -v t="$total_reads" 'BEGIN {printf "%.6f", r/t}')
        else
            frip=0 # Set FRiP to 0 if total reads is 0 or NA
        fi
    elif [[ ! -s "$final_peaks" ]]; then
        log "Warning: Final peak file (${final_peaks}) is empty. Setting FRiP metrics to 0/NA."
        reads_in_peaks=0
        frip=0
    else
        log "Warning: Chr-prefixed BAM ($chr_bam) not found or final peaks empty. Skipping FRiP calculation."
        reads_in_peaks="NA"
        frip="NA"
    fi

    # Calculate peak statistics if final peaks file exists and is not empty
    if [[ -s "$final_peaks" ]]; then
        peak_count=$(wc -l < "$final_peaks" | awk '{print $1}')
        if [[ "$peak_count" -gt 0 ]]; then
            # Use awk to calculate mean length and fold enrichment safely
            read mean_length mean_fold < <(awk '
                BEGIN { sum_len=0; sum_fold=0 }
                { sum_len += ($3 - $2); sum_fold += $7 }
                END {
                    if (NR > 0) {
                        printf "%.2f\t%.2f\n", sum_len/NR, sum_fold/NR
                    } else {
                        print "0.00\t0.00"
                    }
                }' "$final_peaks")
        fi
    else
        log "Warning: Final peak file (${final_peaks}) is empty or does not exist. Setting peak stats to 0."
        peak_count=0
        mean_length=0
        mean_fold=0
    fi


    # Save QC metrics using printf for robustness
    log "Saving QC metrics to ${metrics_csv}..."
    {
        printf "Metric,Value\n"
        printf "Sample,%s\n" "$sample_name"
        printf "Total_Reads,%s\n" "$total_reads"
        printf "Reads_in_Peaks,%s\n" "$reads_in_peaks"
        printf "FRiP,%s\n" "$frip"
        printf "Final_Peak_Count,%s\n" "$peak_count"
        printf "Mean_Peak_Length,%.2f\n" "$mean_length"
        printf "Mean_Fold_Enrichment,%.2f\n" "$mean_fold"
    } > "$metrics_csv" || die "Failed to write QC metrics CSV for $sample_name"
    log "QC metrics calculated and saved for ${sample_name}."

    # Step 9: Generate Peak Statistics File
    log "Step 9: Generating peak statistics file (${stats_txt})..."
    if [[ "$peak_count" -gt 0 ]]; then
        awk 'BEGIN {OFS="\t"; print "Length", "FoldEnrichment"} { print $3-$2, $7 }' "$final_peaks" > "$stats_txt" \
            || die "Failed to generate peak statistics file for $sample_name"
    else
        printf "Length\tFoldEnrichment\n" > "$stats_txt" # Create header only for empty file
    fi
    log "Peak statistics file generated for ${sample_name}."

    # --- Final Summary ---
    log "--- Summary for ${sample_name} ---"
    log "Final Peak Count: $peak_count"
    log "FRiP Score: $frip"
    log "Output Directory: ${output_dir}"
    log "Metrics CSV: ${metrics_csv}"
    log "Final Peaks: ${final_peaks}"
    log "Visualization Peaks: ${viz_peaks}"
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
# bash 5_peak_calling_strict_v2.sh 0

# Parallel using GNU Parallel:
# N=6; seq 0 5 | parallel -j $N bash 5_peak_calling_strict_v2.sh {}

# Parallel using xargs:
# N=6; printf "%s\n" {0..5} | xargs -I {} -P $N bash 5_peak_calling_strict_v2.sh {}