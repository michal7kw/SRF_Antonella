import pandas as pd
import numpy as np
import os

# Define file paths relative to the script location or project root
# Assuming the script is run from d:/Github/SRF_H2AK119Ub_cross_V5
path_pink = '0_check_dysregulation_of_yaf_enriched_genes/YAF_enrichment_results/H2AK119Ub_YAF2_enriched_promoters.csv' # Corrected path
path_orange = '1_find_gene_lists_intersections/output/YAF_SOX.csv'
path_blue = '2_check_dysregulation_of_the_intersected_genes/YAF_dysregulation_analysis_SES/gene_dysregulation_results.csv'

# --- Load Gene Lists ---

# Function to load genes, handling potential errors and stripping versions
def load_gene_set(filepath, gene_col='gene_id', filters=None, remove_version=True):
    """Loads gene IDs from a CSV file into a set, applying optional filters."""
    print(f"Loading {os.path.basename(filepath)}...")
    if not os.path.exists(filepath):
        print(f"Error: File not found at {filepath}")
        return set()
    try:
        df = pd.read_csv(filepath)
        print(f"  Columns: {', '.join(df.columns)}")

        # Find the correct gene ID column (case-insensitive check)
        actual_gene_col = None
        for col in df.columns:
            if col.lower() == gene_col.lower():
                actual_gene_col = col
                break
        if not actual_gene_col:
             # Try common alternatives if primary not found
            common_alternatives = ['geneid', 'gene_name', 'row.names']
            for alt_col_base in common_alternatives:
                 for col in df.columns:
                    if col.lower() == alt_col_base.lower():
                        actual_gene_col = col
                        print(f"  Warning: Using alternative gene column '{actual_gene_col}' for {os.path.basename(filepath)}")
                        break
                 if actual_gene_col: break

        if not actual_gene_col:
            print(f"Error: Gene ID column '{gene_col}' (or alternatives) not found in {os.path.basename(filepath)}")
            return set()
        print(f"  Using gene column: '{actual_gene_col}'")

        # Apply filters if provided
        if filters:
            print(f"  Applying filters: {filters}")
            for col, condition, value in filters:
                 # Check if filter column exists
                if col not in df.columns:
                    print(f"  Warning: Filter column '{col}' not found. Skipping filter.")
                    continue
                try:
                    # Ensure the column is numeric before comparison if possible
                    if pd.api.types.is_numeric_dtype(df[col]):
                        # Convert value to the same type as the column for safe comparison
                        # Handle potential NaN values in the column before getting type
                        first_valid_index = df[col].first_valid_index()
                        if first_valid_index is not None:
                            value_type = type(df[col].loc[first_valid_index])
                            try:
                                comparison_value = value_type(value)
                            except (ValueError, TypeError):
                                 print(f"  Warning: Could not convert filter value '{value}' to type {value_type} for column '{col}'. Skipping filter.")
                                 continue
                        else: # Column is all NaN
                             print(f"  Warning: Column '{col}' contains only NaN values. Skipping filter.")
                             continue

                    else: # Assume string comparison or other types
                        comparison_value = value

                    # Apply the filter safely
                    if condition == '<':
                        df = df[df[col] < comparison_value]
                    elif condition == '>':
                        df = df[df[col] > comparison_value]
                    elif condition == '<=':
                        df = df[df[col] <= comparison_value]
                    elif condition == '>=':
                        df = df[df[col] >= comparison_value]
                    elif condition == '==':
                        df = df[df[col] == comparison_value]
                    else:
                         print(f"  Warning: Unsupported filter condition '{condition}'. Skipping filter.")
                except Exception as e:
                    print(f"  Warning: Could not apply filter '{col} {condition} {value}'. Error: {e}. Skipping filter.")
            print(f"  Genes after filtering: {len(df)}")

        # Extract gene IDs and remove potential version numbers (e.g., .1, .2)
        gene_ids = df[actual_gene_col].dropna().astype(str)
        if remove_version:
             gene_ids = gene_ids.str.replace(r'\.\d+$', '', regex=True)

        gene_set = set(gene_ids)
        print(f"  Loaded {len(gene_set)} unique gene IDs.")
        return gene_set

    except Exception as e:
        print(f"Error loading file {filepath}: {e}")
        return set()

# 1. Pink Circle (H2AK119Ub_SES-YAF2_Promoters)
# User mentioned FDR < 0.01. Let's check if 'padj' or 'qValue' exists and apply filter.
# Assuming 'gene_id' is the column.
# Note: The filename implies enrichment (log2FC > 1) might already be applied.
# We will add the FDR filter based on user feedback.
# Check common FDR columns: padj, qValue, FDR
# pink_filters = [('padj', '<', 0.01)] # Removing filter as 'padj' column not present in this file
print("\n--- Loading Pink Set ---")
# Assuming the file already contains the significantly enriched promoters as per its name
set_pink = load_gene_set(path_pink, gene_col='gene_id', filters=None)

# 2. Orange Circle (SOX2 target genes)
# Assuming 'gene_id' column from YAF_SOX.csv
print("\n--- Loading Orange Set ---")
set_orange = load_gene_set(path_orange, gene_col='gene_id')

# 3. Blue Circle (Genes down by SES_YAF2)
# Filter gene_dysregulation_results.csv for log2FoldChange < 0 and padj < 0.05
# Need to confirm actual column names ('gene_id'?, 'log2FoldChange', 'padj')
blue_filters = [('log2FoldChange', '<', 0), ('padj', '<', 0.05)]
print("\n--- Loading Blue Set ---")
set_blue = load_gene_set(path_blue, gene_col='gene_id', filters=blue_filters) # Corrected gene column to 'gene_id'

# --- Calculate Venn Diagram Segments ---
print("\nCalculating Venn diagram segments...")

# Check if sets were loaded correctly
if not set_pink or not set_orange or not set_blue:
    print("Error: One or more gene sets failed to load. Cannot calculate intersections.")
else:
    # Individual sets (only in that set)
    pink_only = len(set_pink - set_orange - set_blue)
    orange_only = len(set_orange - set_pink - set_blue)
    blue_only = len(set_blue - set_pink - set_orange)

    # Two-way intersections (only in those two)
    pink_orange_only = len((set_pink & set_orange) - set_blue)
    pink_blue_only = len((set_pink & set_blue) - set_orange)
    orange_blue_only = len((set_orange & set_blue) - set_pink)

    # Three-way intersection
    all_three = len(set_pink & set_orange & set_blue)

    # --- Print Results ---
    print("\n--- Calculated Venn Diagram Counts ---")
    print(f"Pink Only (H2AK119Ub_SES-YAF2_Promoters): {pink_only}")
    print(f"Orange Only (SOX2 target genes):          {orange_only}")
    print(f"Blue Only (Genes down by SES_YAF2):       {blue_only}")
    print(f"Pink & Orange Only:                       {pink_orange_only}")
    print(f"Pink & Blue Only:                         {pink_blue_only}")
    print(f"Orange & Blue Only:                       {orange_blue_only}")
    print(f"Pink & Orange & Blue (Center):            {all_three}")
    print("--------------------------------------")

    # Optional: Print total in each set for verification
    print(f"\nTotal in Pink set:   {len(set_pink)}")
    print(f"Total in Orange set: {len(set_orange)}")
    print(f"Total in Blue set:   {len(set_blue)}")

    # Compare the calculated counts above with the Venn diagram image provided.
