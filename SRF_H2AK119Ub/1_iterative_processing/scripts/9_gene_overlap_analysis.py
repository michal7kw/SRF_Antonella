import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np

# File paths
yaf_genes_file = "./analysis/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt"
yaf_full_file = "./analysis/gene_lists_broad/YAF_enriched_genes_broad_full.csv"
sox2_genes_file = "./sox2_binding.csv"
output_dir = "analysis/overlap_analysis"

# create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Read gene lists
def read_gene_list(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def get_regulatory_region_genes(file_path):
    """Extract genes that are enriched in promoter, 5' UTR, and first exon regions"""
    df = pd.read_csv(file_path)
    regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR']
    first_exon = df['annotation'].str.contains('exon 1 of', na=False)
    regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
    return set(df[regulatory_mask]['SYMBOL'].unique())

# Read both gene lists
yaf_genes = read_gene_list(yaf_genes_file)
sox2_genes = read_gene_list(sox2_genes_file)
yaf_regulatory_genes = get_regulatory_region_genes(yaf_full_file)

# Find overlapping and unique genes
overlapping_genes = yaf_genes.intersection(sox2_genes)
yaf_only_genes = yaf_genes.difference(sox2_genes)

# Find overlapping and unique genes for regulatory regions
regulatory_overlapping_genes = yaf_regulatory_genes.intersection(sox2_genes)
regulatory_yaf_only_genes = yaf_regulatory_genes.difference(sox2_genes)

# Print basic statistics
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

# Create Venn diagrams
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

# Save gene lists to files for GO analysis
gene_lists = {
    'all_overlapping_genes.txt': overlapping_genes,
    'all_yaf_only_genes.txt': yaf_only_genes,
    'regulatory_overlapping_genes.txt': regulatory_overlapping_genes,
    'regulatory_yaf_only_genes.txt': regulatory_yaf_only_genes
}

for filename, gene_set in gene_lists.items():
    with open(os.path.join(output_dir, filename), 'w') as f:
        f.write('\n'.join(sorted(gene_set)))

# Calculate enrichment scores for regulatory regions
def calculate_enrichment_scores():
    df = pd.read_csv(yaf_full_file)
    regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR']
    first_exon = df['annotation'].str.contains('exon 1 of', na=False)
    regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
    
    regulatory_df = df[regulatory_mask].copy()
    regulatory_df['is_sox2_target'] = regulatory_df['SYMBOL'].isin(sox2_genes)
    
    # Calculate mean fold change for SOX2 targets vs non-targets
    enrichment_stats = regulatory_df.groupby('is_sox2_target')['fold_change'].agg(['mean', 'std', 'count']).round(3)
    enrichment_stats.to_csv(os.path.join(output_dir, 'regulatory_regions_enrichment_stats.csv'))
    print("\nEnrichment statistics for regulatory regions:")
    print(enrichment_stats)
    
    # Create boxplot of fold changes
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=regulatory_df, x='is_sox2_target', y='fold_change')
    plt.title('H2AK119Ub enrichment in regulatory regions\nSOX2 targets vs non-targets')
    plt.xlabel('Is SOX2 target')
    plt.ylabel('Fold change')
    plt.savefig(os.path.join(output_dir, 'regulatory_regions_enrichment_boxplot.png'))
    plt.close()

calculate_enrichment_scores()

print("\nAnalysis complete. Files generated in", output_dir + ":")
print("1. venn_diagrams.png - Visualization of gene overlaps")
print("2. *_genes.txt - Lists of genes in different categories (ready for GO analysis)")
print("3. regulatory_regions_enrichment_stats.csv - Enrichment statistics for regulatory regions")
print("4. regulatory_regions_enrichment_boxplot.png - Visualization of enrichment scores")
