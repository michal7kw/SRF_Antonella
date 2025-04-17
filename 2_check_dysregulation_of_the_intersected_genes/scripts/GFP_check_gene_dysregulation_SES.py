#!/usr/bin/env python3
"""
Script to check the type of dysregulation of genes listed in GFP_SES.csv 
and GFP_SES_strict.csv by comparing with differential expression results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Configuration ---
# BASE_PATH = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
BASE_PATH = "D:/Github/SRF_H2AK119Ub_cross_V5"

diff_expr_file = F'{BASE_PATH}/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv'

configs = [
    {
        "name": "Non-Strict",
        "gene_list_path": f'{BASE_PATH}/1_find_gene_lists_intersections/output/GFP_SES.csv',
        "output_dir": f'{BASE_PATH}/2_check_dysregulation_of_the_intersected_genes/GFP_dysregulation_analysis_SES'
    },
    {
        "name": "Strict",
        "gene_list_path": f'{BASE_PATH}/1_find_gene_lists_intersections/output/GFP_SES_strict.csv',
        "output_dir": f'{BASE_PATH}/2_check_dysregulation_of_the_intersected_genes/GFP_dysregulation_analysis_SES_strict'
    }
]

# --- Load Differential Expression Data Once ---
print(f"Loading differential expression results from {diff_expr_file}...")
diff_expr_full = pd.read_csv(diff_expr_file, index_col=0)
diff_expr_full['gene_id_base'] = diff_expr_full['gene_id'].str.split('.').str[0]
print(f"Loaded {len(diff_expr_full)} genes from differential expression results.")


# --- Analysis Function ---
def analyze_dysregulation(config_name, gene_list_file, output_dir, diff_expr_data):
    """Performs dysregulation analysis for a given gene list and output directory."""
    
    print(f"--- Starting Analysis for {config_name} ---")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load the gene list
    print(f"Loading gene list from {gene_list_file}...")
    try:
        gene_list = pd.read_csv(gene_list_file)
        print(f"Loaded {len(gene_list)} genes from gene list.")
    except FileNotFoundError:
        print(f"Error: Gene list file not found at {gene_list_file}. Skipping this configuration.")
        return

    # Extract gene_id without version
    gene_list['gene_id_base'] = gene_list['gene_id'].str.split('.').str[0]

    # Merge the datasets
    print("Matching genes between the two datasets...")
    # Make a copy to avoid modifying the original diff_expr_data DataFrame
    diff_expr = diff_expr_data.copy() 
    merged_data = pd.merge(gene_list, diff_expr, 
                          left_on='gene_id_base', right_on='gene_id_base', 
                          how='left', suffixes=('', '_diff'))

    # Define dysregulation categories
    # Using standard thresholds: |log2FC| > 1 and padj < 0.05
    def categorize_dysregulation(row):
        if pd.isna(row['log2FoldChange']) or pd.isna(row['padj']):
            return "Not in differential expression data"
        elif row['padj'] >= 0.05:
            return "Not significantly dysregulated"
        elif row['log2FoldChange'] > 1:
            return "Up-regulated"
        elif row['log2FoldChange'] < -1:
            return "Down-regulated"
        else:
            return "Minimally changed"

    # Apply categorization
    merged_data['dysregulation'] = merged_data.apply(categorize_dysregulation, axis=1)

    # Count genes in each category
    dysreg_counts = merged_data['dysregulation'].value_counts()
    print(f"Dysregulation summary ({config_name}):")
    print(dysreg_counts)

    # Save the results
    output_file = os.path.join(output_dir, 'gene_dysregulation_results.csv')
    merged_data.to_csv(output_file, index=False) # Avoid writing pandas index
    print(f"Detailed results saved to {output_file}")

    # Create a bar plot of dysregulation categories
    plt.figure(figsize=(10, 6))
    sns.barplot(x=dysreg_counts.index, y=dysreg_counts.values)
    plt.title(f'Dysregulation Categories for GFP_SES ({config_name}) Genes')
    plt.xlabel('Dysregulation Category')
    plt.ylabel('Number of Genes')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    bar_plot_path = os.path.join(output_dir, 'dysregulation_categories.png')
    plt.savefig(bar_plot_path)
    plt.close() # Close the plot to free memory

    # Create a volcano plot
    plt.figure(figsize=(10, 8))
    # Filter to only genes in our list that have differential expression data
    volcano_data = merged_data.dropna(subset=['log2FoldChange', 'padj']).copy() # Use .copy() to avoid SettingWithCopyWarning

    # Create a new column for point colors
    volcano_data['color'] = 'gray'  # Default color
    volcano_data.loc[(volcano_data['padj'] < 0.05) & (volcano_data['log2FoldChange'] > 1), 'color'] = 'red'
    volcano_data.loc[(volcano_data['padj'] < 0.05) & (volcano_data['log2FoldChange'] < -1), 'color'] = 'blue'

    # Plot
    plt.scatter(
        volcano_data['log2FoldChange'], 
        -np.log10(volcano_data['padj']), 
        c=volcano_data['color'],
        alpha=0.6
    )

    # Add labels for top dysregulated genes
    top_genes = volcano_data.sort_values('padj').head(20)
    for _, row in top_genes.iterrows():
        # Check if gene_name exists and is not NaN before plotting text
        if 'gene_name' in row and pd.notna(row['gene_name']):
            plt.text(row['log2FoldChange'], -np.log10(row['padj']), row['gene_name'], 
                     fontsize=8, ha='center', va='bottom')

    # Add lines for thresholds
    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.3)
    plt.axvline(1, color='gray', linestyle='--', alpha=0.3)
    plt.axvline(-1, color='gray', linestyle='--', alpha=0.3)

    plt.title(f'Volcano Plot of GFP_SES ({config_name}) Genes')
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(adjusted p-value)')
    plt.tight_layout()
    volcano_plot_path = os.path.join(output_dir, 'volcano_plot.png')
    plt.savefig(volcano_plot_path)
    plt.close() # Close the plot

    # Create a heatmap of the top dysregulated genes
    # Get top 30 most significantly dysregulated genes that have a gene name
    top_dysreg = merged_data.dropna(subset=['padj', 'gene_name']).sort_values('padj').head(30)
    
    if not top_dysreg.empty:
        plt.figure(figsize=(12, max(10, len(top_dysreg) * 0.4))) # Adjust height based on number of genes
        heatmap_data = top_dysreg[['gene_name', 'log2FoldChange', 'padj']].set_index('gene_name')
        sns.heatmap(heatmap_data[['log2FoldChange', 'padj']], cmap='coolwarm', annot=True, fmt='.2g') # Select only numeric columns for heatmap
        plt.title(f'Top {len(top_dysreg)} Significantly Dysregulated Genes ({config_name})')
        plt.tight_layout()
        heatmap_plot_path = os.path.join(output_dir, 'top_dysregulated_heatmap.png')
        plt.savefig(heatmap_plot_path)
        plt.close() # Close the plot
        heatmap_generated = True
    else:
        print("No significantly dysregulated genes with names found to generate heatmap.")
        heatmap_generated = False
        heatmap_plot_path = "N/A"


    print(f"Analysis complete for {config_name}. Generated the following visualizations:")
    print(f"1. {bar_plot_path}")
    print(f"2. {volcano_plot_path}")
    if heatmap_generated:
        print(f"3. {heatmap_plot_path}")
    else:
        print("3. Top dysregulated heatmap: Not generated (no suitable data).")


# --- Main Loop ---
for config in configs:
    analyze_dysregulation(
        config_name=config["name"],
        gene_list_file=config["gene_list_path"],
        output_dir=config["output_dir"],
        diff_expr_data=diff_expr_full
    )

print("--- All analyses finished ---")
