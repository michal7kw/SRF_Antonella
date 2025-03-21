#!/usr/bin/env python3

import sys
import re

def convert_scientific_to_decimal(match):
    """Convert scientific notation to decimal format."""
    number = float(match.group(0))
    return str(int(number))

def process_file(input_file, output_file):
    """Process the broadPeak file and convert scientific notation."""
    
    # Regular expression to match scientific notation
    # Matches numbers like 7.5e+07, 1.23E-04, etc.
    scientific_pattern = r'\d+\.?\d*[eE][+-]?\d+'
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # Convert any scientific notation numbers in the line
                corrected_line = re.sub(scientific_pattern, convert_scientific_to_decimal, line)
                outfile.write(corrected_line)
                
    except FileNotFoundError:
        print(f"Error: Could not find file {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python fix_scientific_notation_v2.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = input_file.rsplit('.', 1)[0] + '_fixed.' + input_file.rsplit('.', 1)[1]
    
    process_file(input_file, output_file)
    print(f"Processed file saved as: {output_file}")

if __name__ == "__main__":
    main() 