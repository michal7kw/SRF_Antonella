#!/usr/bin/env python3
"""
Script to merge multiple broadPeak files from Cut&Tag experiments.
"""
import os
import sys
import glob

def merge_broadpeak_files(input_files, output_file):
    """
    Merge multiple broadPeak files into a single file.
    
    Args:
        input_files (list): List of input broadPeak files
        output_file (str): Path to output merged broadPeak file
    """
    print(f"Merging {len(input_files)} files into {output_file}")
    
    # Check if all input files exist
    for file in input_files:
        if not os.path.exists(file):
            print(f"Error: Input file {file} does not exist")
            sys.exit(1)
    
    # Open output file for writing
    with open(output_file, 'w') as outf:
        # Process each input file
        for file in input_files:
            print(f"Processing {file}")
            try:
                with open(file, 'r') as inf:
                    # Copy all lines from input file to output file
                    for line in inf:
                        outf.write(line)
            except Exception as e:
                print(f"Error processing {file}: {e}")
                sys.exit(1)
    
    print(f"Successfully merged files into {output_file}")

def main():
    # Directory containing the broadPeak files
    directory = os.path.dirname(os.path.abspath(__file__))
    
    # Pattern for GFP broadPeak files
    pattern = os.path.join(directory, "GFP_*_broad_peaks_final.broadPeak")
    
    # Find all matching files
    input_files = sorted(glob.glob(pattern))
    
    if not input_files:
        print("Error: No GFP broadPeak files found")
        sys.exit(1)
    
    print(f"Found {len(input_files)} files: {', '.join(input_files)}")
    
    # Output file
    output_file = os.path.join(directory, "GFP.broadPeak")
    
    # Merge the files
    merge_broadpeak_files(input_files, output_file)

if __name__ == "__main__":
    main()
