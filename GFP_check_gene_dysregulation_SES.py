#!/usr/bin/env python3
"""
Script to check the type of dysregulation of genes listed in YAF_SES.csv
by comparing with differential expression results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# File paths
gene_list_file = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/YAF_SES_strict.csv'
output_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/dysregulation_analysis_SES_strict'

# gene_list_file = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/YAF_SES.csv'
# output_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/dysregulation_analysis_SES'

diff_expr_file = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv'



# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the gene list
print(f"Loading gene list from {gene_list_file}...")
gene_list = pd.read_csv(gene_list_file)
print(f"Loaded {len(gene_list)} genes from gene list.")

# Load the differential expression results
print(f"Loading differential expression results from {diff_expr_file}...")
diff_expr = pd.read_csv(diff_expr_file, index_col=0)
print(f"Loaded {len(diff_expr)} genes from differential expression results.")

# Extract gene_id without version from both datasets
gene_list['gene_id_base'] = gene_list['gene_id'].str.split('.').str[0]
diff_expr['gene_id_base'] = diff_expr['gene_id'].str.split('.').str[0]

# Merge the datasets
print("Matching genes between the two datasets...")
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
print("\nDysregulation summary:")
print(dysreg_counts)

# Save the results
output_file = os.path.join(output_dir, 'gene_dysregulation_results.csv')
merged_data.to_csv(output_file)
print(f"\nDetailed results saved to {output_file}")

# Create a bar plot of dysregulation categories
plt.figure(figsize=(10, 6))
sns.barplot(x=dysreg_counts.index, y=dysreg_counts.values)
plt.title('Dysregulation Categories for YAF_SES Genes')
plt.xlabel('Dysregulation Category')
plt.ylabel('Number of Genes')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'dysregulation_categories.png'))

# Create a volcano plot
plt.figure(figsize=(10, 8))
# Filter to only genes in our list that have differential expression data
volcano_data = merged_data.dropna(subset=['log2FoldChange', 'padj'])

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
    plt.text(row['log2FoldChange'], -np.log10(row['padj']), row['gene_name'], 
             fontsize=8, ha='center', va='bottom')

# Add lines for thresholds
plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.3)
plt.axvline(1, color='gray', linestyle='--', alpha=0.3)
plt.axvline(-1, color='gray', linestyle='--', alpha=0.3)

plt.title('Volcano Plot of YAF_SES Genes')
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(adjusted p-value)')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'volcano_plot.png'))

# Create a heatmap of the top dysregulated genes
# Get top 30 most significantly dysregulated genes
top_dysreg = merged_data.dropna(subset=['padj']).sort_values('padj').head(30)
plt.figure(figsize=(12, 10))
heatmap_data = top_dysreg[['gene_name', 'log2FoldChange', 'padj']].set_index('gene_name')
sns.heatmap(heatmap_data, cmap='coolwarm', annot=True, fmt='.2g')
plt.title('Top 30 Dysregulated Genes')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'top_dysregulated_heatmap.png'))

print("\nAnalysis complete. Generated the following visualizations:")
print(f"1. {os.path.join(output_dir, 'dysregulation_categories.png')}")
print(f"2. {os.path.join(output_dir, 'volcano_plot.png')}")
print(f"3. {os.path.join(output_dir, 'top_dysregulated_heatmap.png')}")
