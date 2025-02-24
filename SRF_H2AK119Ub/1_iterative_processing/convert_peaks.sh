#!/bin/bash
#SBATCH --job-name=convert_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/convert_peaks.err"
#SBATCH --output="logs/convert_peaks.out"

# Description:
# This script converts chromosome names in broadPeak files by removing the "chr" prefix.
# This is necessary for compatibility with certain downstream analysis tools.
# It also creates a backup of the original peak files before conversion.

# Input:
# - broadPeak files located in the ./analysis/5_peak_calling directory, named *_broad_peaks_final.broadPeak

# Output:
# - Modified broadPeak files in the same ./analysis/5_peak_calling directory, with "chr" prefix removed from chromosome names.
# - Backup of original broadPeak files in ./analysis/5_peak_calling/original_peaks_backup directory.

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Directory containing peak files
PEAKS_DIR="./analysis/5_peak_calling"

# Create a backup directory
BACKUP_DIR="${PEAKS_DIR}/original_peaks_backup"
mkdir -p "${BACKUP_DIR}"

# Process each broadPeak file
for peakfile in "${PEAKS_DIR}"/*_broad_peaks_final.broadPeak; do
    # Get just the filename
    filename=$(basename "${peakfile}")
    
    # Create backup
    cp "${peakfile}" "${BACKUP_DIR}/${filename}"
    
    # Convert chromosome names and save to temporary file
    sed 's/^chr//' "${peakfile}" > "${peakfile}.tmp"
    
    # Move temporary file to original location
    mv "${peakfile}.tmp" "${peakfile}"
    
    echo "Converted ${filename}"
done

echo "All peak files have been converted and originals backed up in ${BACKUP_DIR}"