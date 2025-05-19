# This script performs overlap analysis between YAF-enriched genes and SOX2 target genes.

#### Input files: ####
# - analysis/gene_lists_broad/YAF_enriched_genes
#     List of genes enriched in YAF samples (gene symbols)
# ---------------------------------------------------------------------------------------------------------
# MIR8071-2
# MIR8071-1
# TBX4
# TBC1D3P1-DHX40P1
# ATP9B
# MST1L
# ---------------------------------------------------------------------------------------------------------

# - analysis/gene_lists_broad/YAF_enriched_genes
#     Full data for YAF-enriched genes including annotation and fold change
# ---------------------------------------------------------------------------------------------------------
# "ENTREZID","SYMBOL","distanceToTSS","annotation","fold_change"
# "102466889","MIR8071-2",6171,"Exon (ENST00000497397.1/ENST00000497397.1, exon 2 of 3)",1.76251149051725
# "102465871","MIR8071-1",-31123,"Exon (ENST00000497872.4/ENST00000497872.4, exon 1 of 5)",1.63204716725618
# "9496","TBX4",-7081,"Distal Intergenic",1.6301302738487
# ---------------------------------------------------------------------------------------------------------

# - sox2_binding.csv:
#     List of SOX2 target genes
# ---------------------------------------------------------------------------------------------------------
# A2M
# A4GALT
# AADAC
# AADAT
# AAK1
# AAMDC
# AANAT
# ---------------------------------------------------------------------------------------------------------

#### Output files (in analysis/11_gene_overlap_analysis/): ####
# 1. venn_diagrams.png:
#     - Two Venn diagrams showing overlap between YAF and SOX2 genes
#     - One for all genes, one for genes in regulatory regions only
# 2. Gene lists for GO analysis:
#     - all_overlapping_genes.txt: Genes found in both YAF and SOX2 sets
#     - all_yaf_only_genes.txt: Genes unique to YAF set
#     - regulatory_overlapping_genes.txt: Genes in regulatory regions found in both sets
#     - regulatory_yaf_only_genes.txt: Genes in regulatory regions unique to YAF
# 3. regulatory_regions_enrichment_stats.csv:
#     Statistics comparing fold changes between SOX2 targets and non-targets
# 4. regulatory_regions_enrichment_boxplot.png:
#     Visualization of enrichment scores in regulatory regions

import pandas as pd
import os
import matplotlib
matplotlib.use('Agg') # Add this line to use a non-interactive backend
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np
import sys

# Define file paths and create output directory
OUTPUT_DIR = sys.argv[1]
# BASE_PATH = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
# BASE_PATH = "D:/Github/SRF_H2AK119Ub_cross_V5"
BASE_PATH = "/mnt/d/Github/SRF_H2AK119Ub_cross_V5"

YAF_ALL_ENRICHED_GENES_FILE =  os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_significant/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt")
YAF_PROMOTER_GENES_FILE = os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_significant/gene_lists_broad/YAF_enriched_genes_broad_promoters.txt")
YAF_FULL_FILE =  os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_significant/gene_lists_broad/YAF_enriched_genes_broad_full.csv") # Still needed for calculate_enrichment_scores
SOX2_GENES_FILE =  os.path.join(BASE_PATH, "COMMON_DATA/sox2_binding.csv")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def read_gene_list(file_path):
    """Read a list of genes from a file, one gene per line"""
    try:
        with open(file_path, 'r', encoding='latin1') as f:
            return set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def get_regulatory_region_genes(file_path):
    """
    Extract genes that are enriched in regulatory regions for enrichment score calculation:
    - Promoter regions (up to 3kb upstream)
    - 5' UTR
    - First exon
    
    Returns a set of gene symbols. This function is used by calculate_enrichment_scores.
    """
    try:
        df = pd.read_csv(file_path)
        regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR', 'Exon']
        first_exon = df['annotation'].str.contains('exon 1 of', na=False)
        regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
        return set(df[regulatory_mask]['SYMBOL'].unique())
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

# Read gene lists
yaf_all_enriched_genes = read_gene_list(YAF_ALL_ENRICHED_GENES_FILE)
yaf_promoter_genes = read_gene_list(YAF_PROMOTER_GENES_FILE)
sox2_genes = read_gene_list(SOX2_GENES_FILE)

# --- Venn Diagram 1: YAF Promoter Genes vs. SOX2 Target Genes ---
promoter_overlapping_genes = yaf_promoter_genes.intersection(sox2_genes)
promoter_yaf_only_genes = yaf_promoter_genes.difference(sox2_genes)
promoter_sox2_only_genes = sox2_genes.difference(yaf_promoter_genes) # For completeness if needed

# --- Venn Diagram 2: All YAF-Enriched Genes vs. SOX2 Target Genes ---
all_enriched_overlapping_genes = yaf_all_enriched_genes.intersection(sox2_genes)
all_enriched_yaf_only_genes = yaf_all_enriched_genes.difference(sox2_genes)
all_enriched_sox2_only_genes = sox2_genes.difference(yaf_all_enriched_genes) # For completeness

# Print summary statistics
print(f"\nBasic Statistics:")
print(f"Total All YAF-Enriched genes: {len(yaf_all_enriched_genes)}")
print(f"Total YAF Promoter genes: {len(yaf_promoter_genes)}")
print(f"Total SOX2 target genes: {len(sox2_genes)}")

print(f"\nOverlap 1: YAF Promoter Genes vs. SOX2 Target Genes")
print(f"  Number of overlapping genes (YAF Promoter & SOX2): {len(promoter_overlapping_genes)}")
print(f"  Number of YAF Promoter only genes: {len(promoter_yaf_only_genes)}")
print(f"  Number of SOX2 only genes (vs YAF Promoter): {len(promoter_sox2_only_genes)}")

print(f"\nOverlap 2: All YAF-Enriched Genes vs. SOX2 Target Genes")
print(f"  Number of overlapping genes (All YAF-Enriched & SOX2): {len(all_enriched_overlapping_genes)}")
print(f"  Number of All YAF-Enriched only genes: {len(all_enriched_yaf_only_genes)}")
print(f"  Number of SOX2 only genes (vs All YAF-Enriched): {len(all_enriched_sox2_only_genes)}")

# Create and save Venn diagrams
plt.figure(figsize=(16, 8)) # Adjusted figure size

# Diagram 1: YAF Promoter Binding vs. SOX2 Target Genes
plt.subplot(1, 2, 1)
venn2_diagram1 = venn2(
    (len(promoter_yaf_only_genes), len(promoter_sox2_only_genes), len(promoter_overlapping_genes)),
    set_labels=('YAF Promoter Genes', 'SOX2 Target Genes'),
    set_colors=('#F9A252', '#64A0C8'), alpha=0.7
)
plt.title('YAF Promoter Genes vs. SOX2 Target Genes', fontsize=14)

# Diagram 2: All YAF-Enriched Genes vs. SOX2 Target Genes
plt.subplot(1, 2, 2)
venn2_diagram2 = venn2(
    (len(all_enriched_yaf_only_genes), len(all_enriched_sox2_only_genes), len(all_enriched_overlapping_genes)),
    set_labels=('All YAF-Enriched Genes', 'SOX2 Target Genes'),
    set_colors=('#F9A252', '#64A0C8'), alpha=0.7
)
plt.title('All YAF-Enriched Genes vs. SOX2 Target Genes', fontsize=14)

plt.tight_layout(rect=(0, 0, 1, 0.96)) # Add rect to make space for suptitle
plt.suptitle("YAF and SOX2 Gene Overlap Analysis", fontsize=16, y=0.99)
plt.savefig(os.path.join(OUTPUT_DIR, 'venn_diagrams.png'))
plt.close()

# Save gene lists for downstream GO analysis
gene_lists = {
    'promoter_overlapping_genes.txt': promoter_overlapping_genes,
    'promoter_yaf_only_genes.txt': promoter_yaf_only_genes,
    'all_overlapping_genes.txt': all_enriched_overlapping_genes, # Corresponds to "all_overlapping_genes.txt" in shell script
    'all_yaf_only_genes.txt': all_enriched_yaf_only_genes      # Corresponds to "all_yaf_only_genes.txt" in shell script
}
# Note: The shell script also mentions regulatory_overlapping_genes.txt and regulatory_yaf_only_genes.txt.
# These are now promoter_overlapping_genes.txt and promoter_yaf_only_genes.txt.

for filename, gene_set in gene_lists.items():
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        with open(filepath, 'w') as f:
            f.write('\n'.join(sorted(gene_set)))
    except Exception as e:
        print(f"Error writing to file {filepath}: {e}")
        sys.exit(1)

def calculate_enrichment_scores():
    """
    Calculate and visualize enrichment scores for regulatory regions.
    Compares H2AK119Ub fold changes between SOX2 targets and non-targets.
    
    Outputs:
    - CSV file with mean, std, and count statistics
    - Two boxplots visualization of fold changes (with and without outliers)
    """
    try:
        df = pd.read_csv(YAF_FULL_FILE)
        regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR', 'Exon']
        first_exon = df['annotation'].str.contains('exon 1 of', na=False)
        regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
        
        regulatory_df = df[regulatory_mask].copy()
        regulatory_df['is_sox2_target'] = regulatory_df['SYMBOL'].isin(sox2_genes)
        
        # Calculate statistics
        enrichment_stats = regulatory_df.groupby('is_sox2_target')['fold_change'].agg(['mean', 'std', 'count']).round(3)
        enrichment_stats_path = os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_stats.csv')
        enrichment_stats.to_csv(enrichment_stats_path)
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
        enrichment_boxplot_path = os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_boxplot.png')
        plt.savefig(enrichment_boxplot_path, bbox_inches='tight')
        plt.close()
    
    except FileNotFoundError:
        print(f"Error: The file {YAF_FULL_FILE} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during enrichment score calculation: {e}")
        sys.exit(1)

calculate_enrichment_scores()

# Print summary of generated files
print("\nAnalysis complete. Files generated in", OUTPUT_DIR + ":")
print("1. venn_diagrams.png - Visualization of gene overlaps")
print("2. *_genes.txt - Lists of genes in different categories (ready for GO analysis)")
print("3. regulatory_regions_enrichment_stats.csv - Enrichment statistics for regulatory regions")
print("4. regulatory_regions_enrichment_boxplot.png - Visualization of enrichment scores")
