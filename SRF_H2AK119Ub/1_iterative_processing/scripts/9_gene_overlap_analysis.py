"""
This script performs overlap analysis between YAF-enriched genes and SOX2 target genes.

Input files:
- ./analysis/gene_lists_{peak_type}/YAF_enriched_genes_{peak_type}_symbols.txt:
    List of genes enriched in YAF samples (gene symbols)
- ./analysis/gene_lists_{peak_type}/YAF_enriched_genes_{peak_type}_full.csv:
    Full data for YAF-enriched genes including annotation and fold change
- ./sox2_binding.csv:
    List of SOX2 target genes

Output files (in analysis/overlap_analysis_{peak_type}/):
1. venn_diagrams.png:
    - Two Venn diagrams showing overlap between YAF and SOX2 genes
    - One for all genes, one for genes in regulatory regions only
2. Gene lists for GO analysis:
    - all_overlapping_genes.txt: Genes found in both YAF and SOX2 sets
    - all_yaf_only_genes.txt: Genes unique to YAF set
    - regulatory_overlapping_genes.txt: Genes in regulatory regions found in both sets
    - regulatory_yaf_only_genes.txt: Genes in regulatory regions unique to YAF
3. regulatory_regions_enrichment_stats.csv:
    Statistics comparing fold changes between SOX2 targets and non-targets
4. regulatory_regions_enrichment_boxplot.png:
    Visualization of enrichment scores in regulatory regions
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np

# Define file paths and create output directory
peak_type = "broad"
yaf_genes_file = f"./analysis/gene_lists_{peak_type}/YAF_enriched_genes_{peak_type}_symbols.txt"
yaf_full_file = f"./analysis/gene_lists_{peak_type}/YAF_enriched_genes_{peak_type}_full.csv"
sox2_genes_file = "./sox2_binding.csv"
output_dir = f"analysis/overlap_analysis_{peak_type}"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def read_gene_list(file_path):
    """Read a list of genes from a file, one gene per line"""
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def get_regulatory_region_genes(file_path):
    """
    Extract genes that are enriched in regulatory regions:
    - Promoter regions (up to 3kb upstream)
    - 5' UTR
    - First exon
    
    Returns a set of gene symbols
    """
    df = pd.read_csv(file_path)
    regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR']
    first_exon = df['annotation'].str.contains('exon 1 of', na=False)
    regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
    return set(df[regulatory_mask]['SYMBOL'].unique())

# Read gene lists and find overlaps
yaf_genes = read_gene_list(yaf_genes_file)
sox2_genes = read_gene_list(sox2_genes_file)
yaf_regulatory_genes = get_regulatory_region_genes(yaf_full_file)

# Find overlapping and unique genes for all regions
overlapping_genes = yaf_genes.intersection(sox2_genes)
yaf_only_genes = yaf_genes.difference(sox2_genes)

# Find overlapping and unique genes specifically in regulatory regions
regulatory_overlapping_genes = yaf_regulatory_genes.intersection(sox2_genes)
regulatory_yaf_only_genes = yaf_regulatory_genes.difference(sox2_genes)

# Print summary statistics
print(f"\nBasic Statistics:")
print(f"Total YAF enriched genes: {len(yaf_genes)}")
print(f"YAF genes in regulatory regions: {len(yaf_regulatory_genes)}")
print(f"Total SOX2 target genes: {len(sox2_genes)}")
print(f"\nAll genes overlap:")
print(f"Number of overlapping genes: {len(overlapping_genes)}")
print(f"Number of YAF-only genes: {len(yaf_only_genes)}")
print(f"\nRegulatory regions overlap:")
print(f"Number of overlapping regulatory genes: {len(regulatory_overlapping_genes)}")
print(f"Number of YAF-only regulatory genes: {len(regulatory_yaf_only_genes)}")

# Create and save Venn diagrams
plt.figure(figsize=(15, 7))

plt.subplot(1, 2, 1)
venn2([yaf_genes, sox2_genes], ('YAF enriched genes', 'SOX2 target genes'))
plt.title('All genes overlap')

plt.subplot(1, 2, 2)
venn2([yaf_regulatory_genes, sox2_genes], ('YAF regulatory genes', 'SOX2 target genes'))
plt.title('Regulatory regions overlap')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'venn_diagrams.png'))
plt.close()

# Save gene lists for downstream GO analysis
gene_lists = {
    'all_overlapping_genes.txt': overlapping_genes,
    'all_yaf_only_genes.txt': yaf_only_genes,
    'regulatory_overlapping_genes.txt': regulatory_overlapping_genes,
    'regulatory_yaf_only_genes.txt': regulatory_yaf_only_genes
}

for filename, gene_set in gene_lists.items():
    with open(os.path.join(output_dir, filename), 'w') as f:
        f.write('\n'.join(sorted(gene_set)))

def calculate_enrichment_scores():
    """
    Calculate and visualize enrichment scores for regulatory regions.
    Compares H2AK119Ub fold changes between SOX2 targets and non-targets.
    
    Outputs:
    - CSV file with mean, std, and count statistics
    - Two boxplots visualization of fold changes (with and without outliers)
    """
    df = pd.read_csv(yaf_full_file)
    regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR']
    first_exon = df['annotation'].str.contains('exon 1 of', na=False)
    regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
    
    regulatory_df = df[regulatory_mask].copy()
    regulatory_df['is_sox2_target'] = regulatory_df['SYMBOL'].isin(sox2_genes)
    
    # Calculate statistics
    enrichment_stats = regulatory_df.groupby('is_sox2_target')['fold_change'].agg(['mean', 'std', 'count']).round(3)
    enrichment_stats.to_csv(os.path.join(output_dir, 'regulatory_regions_enrichment_stats.csv'))
    print("\nEnrichment statistics for regulatory regions:")
    print(enrichment_stats)
    
    # Create visualization with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # First plot - with outliers
    sns.boxplot(data=regulatory_df, x='is_sox2_target', y='fold_change', ax=ax1)
    ax1.set_title('With outliers')
    ax1.set_xlabel('Is SOX2 target')
    ax1.set_ylabel('Fold change')
    
    # Second plot - without outliers
    sns.boxplot(data=regulatory_df, x='is_sox2_target', y='fold_change', 
                showfliers=False, ax=ax2)
    ax2.set_title('Without outliers')
    ax2.set_xlabel('Is SOX2 target')
    ax2.set_ylabel('Fold change')
    
    # Add overall title
    plt.suptitle('H2AK119Ub enrichment in regulatory regions\nSOX2 targets vs non-targets', 
                 y=1.05)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulatory_regions_enrichment_boxplot.png'), 
                bbox_inches='tight')
    plt.close()

calculate_enrichment_scores()

# Print summary of generated files
print("\nAnalysis complete. Files generated in", output_dir + ":")
print("1. venn_diagrams.png - Visualization of gene overlaps")
print("2. *_genes.txt - Lists of genes in different categories (ready for GO analysis)")
print("3. regulatory_regions_enrichment_stats.csv - Enrichment statistics for regulatory regions")
print("4. regulatory_regions_enrichment_boxplot.png - Visualization of enrichment scores")
