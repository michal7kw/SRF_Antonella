#!/usr/bin/env python3
"""
Script to merge multiple broadPeak files from Cut&Tag experiments
for specified sample types (e.g., GFP, YAF) within a given directory.
"""
import os
import sys
import glob
import argparse

def merge_broadpeak_files(input_files, output_file):
    """
    Merge multiple broadPeak files into a single file.

    Args:
        input_files (list): List of input broadPeak files
        output_file (str): Path to output merged broadPeak file
    """
    if not input_files:
        print(f"Warning: No input files found for pattern leading to {output_file}. Skipping.")
        return False

    print(f"Merging {len(input_files)} files into {output_file}")

    # Check if all input files exist (glob should ensure this, but good practice)
    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Error: Input file {file_path} does not exist though it was globbed.")
            # Continue, as glob might have found some valid files
            return False # Or handle more gracefully

    try:
        with open(output_file, 'w') as outf:
            for file_path in input_files:
                print(f"Processing {file_path}")
                try:
                    with open(file_path, 'r') as inf:
                        for line in inf:
                            outf.write(line)
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
                    # Optionally, decide if one failed file should stop the whole merge
                    # For now, we'll let it try other files if any, but the output might be incomplete.
                    # Consider sys.exit(1) if one failure is critical.
                    return False # Indicate failure for this merge operation
        print(f"Successfully merged files into {output_file}")
        return True
    except Exception as e:
        print(f"Error opening or writing to output file {output_file}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Merge broadPeak files for GFP and YAF samples from a specified directory.")
    parser.add_argument("peak_directory", help="Directory containing the broadPeak files (e.g., GFP_*_broad_peaks_final.broadPeak). Merged files will also be saved here.")
    
    args = parser.parse_args()
    
    peak_dir = os.path.abspath(args.peak_directory)

    if not os.path.isdir(peak_dir):
        print(f"Error: Provided directory '{peak_dir}' does not exist or is not a directory.")
        sys.exit(1)

    print(f"Target peak directory: {peak_dir}")

    sample_types = ["GFP", "YAF"]
    all_successful = True

    for sample_type in sample_types:
        print(f"\nProcessing sample type: {sample_type}")
        
        # Pattern for broadPeak files
        pattern = os.path.join(peak_dir, f"{sample_type}_*_broad_peaks_final.broadPeak")
        
        # Find all matching files
        input_files = sorted(glob.glob(pattern))
        
        if not input_files:
            print(f"Warning: No {sample_type} broadPeak files found in '{peak_dir}' matching pattern '{os.path.basename(pattern)}'.")
            # This is not necessarily an error for the script's run, just for this sample type.
            # We'll mark all_successful as False if we expect files for every type.
            # For now, let's assume it's okay if one type has no files.
            continue 
        
        print(f"Found {len(input_files)} {sample_type} files: {', '.join(map(os.path.basename, input_files))}")
        
        # Output file
        output_file = os.path.join(peak_dir, f"{sample_type}.broadPeak")
        
        # Merge the files
        if not merge_broadpeak_files(input_files, output_file):
            print(f"Failed to merge files for {sample_type}.")
            all_successful = False # Mark that at least one merge failed

    if all_successful:
        print("\nAll sample types processed successfully.")
    else:
        print("\nSome sample types could not be processed or merged successfully. Please check logs.")
        sys.exit(1) # Exit with error if any merge failed

if __name__ == "__main__":
    main()