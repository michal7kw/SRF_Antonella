# This script performs overlap analysis between YAF-enriched genes, SOX2 target genes,
# and YAF downregulated genes.

#### Input files: ####
# - analysis/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt
#     List of genes enriched in YAF samples (gene symbols) - H2AK119Ub-YAF2 enriched peaks (gene bodies)
# - analysis/gene_lists_broad/YAF_enriched_genes_broad_promoters.txt
#     List of genes enriched in YAF samples (promoters) - H2AK119Ub-YAF2 enriched peaks (promoters)
# - analysis/gene_lists_broad/YAF_enriched_genes_broad_full.csv
#     Full data for YAF-enriched genes including annotation and fold change
# - COMMON_DATA/sox2_binding.csv:
#     List of SOX2 target genes
# - SRF_RNA/results/deseq2/YAF_vs_GFP/summary_files/YAF_vs_GFP/down_regulated.csv
#     List of YAF downregulated genes

#### Output files (in analysis/11_gene_overlap_analysis_extended/): ####
# 1. venn_diagrams_2way.png:
#     - Two Venn diagrams showing overlap between YAF and SOX2 genes (original)
# 2. venn_diagrams_3way.png:
#     - Two Venn diagrams showing 3-way overlap:
#       - H2AK119Ub-YAF2 (gene bodies) vs SOX2 vs YAF downregulated
#       - H2AK119Ub-YAF2 (promoters) vs SOX2 vs YAF downregulated
# 3. Gene lists for GO analysis (original and new 3-way overlaps)
# 4. regulatory_regions_enrichment_stats.csv: (original)
#     Statistics comparing fold changes between SOX2 targets and non-targets
# 5. regulatory_regions_enrichment_boxplot.png: (original)
#     Visualization of enrichment scores in regulatory regions

import pandas as pd
import os
import matplotlib
matplotlib.use('Agg') # Add this line to use a non-interactive backend
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3 # Modified import
import seaborn as sns
import numpy as np
import sys

# Define file paths and create output directory
OUTPUT_DIR = sys.argv[1]
# BASE_PATH = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
# BASE_PATH = "D:/Github/SRF_H2AK119Ub_cross_V5"
BASE_PATH = "/mnt/d/Github/SRF_H2AK119Ub_cross_V5" # Ensure this is correct for your environment

YAF_ALL_ENRICHED_GENES_FILE =  os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_strict_significant/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt")
YAF_PROMOTER_GENES_FILE = os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_strict_significant/gene_lists_broad/YAF_enriched_genes_broad_promoters.txt")
YAF_FULL_FILE =  os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_strict_significant/gene_lists_broad/YAF_enriched_genes_broad_full.csv") # Still needed for calculate_enrichment_scores
SOX2_GENES_FILE =  os.path.join(BASE_PATH, "COMMON_DATA/sox2_binding.csv")
YAF_DOWNREGULATED_GENES_FILE = os.path.join(BASE_PATH, "SRF_RNA/results/deseq2/YAF_vs_GFP/summary_files/YAF_vs_GFP/down_regulated.csv") # New file

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

def read_gene_list_from_csv(file_path, symbol_column='gene_symbol'):
    """Read a list of gene symbols from a CSV file."""
    try:
        df = pd.read_csv(file_path)
        if symbol_column not in df.columns:
            print(f"Error: Column '{symbol_column}' not found in {file_path}.")
            sys.exit(1)
        # Drop NA values in the symbol column to avoid issues with set operations
        return set(df[symbol_column].dropna().unique())
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV file {file_path}: {e}")
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
yaf_downregulated_genes = read_gene_list_from_csv(YAF_DOWNREGULATED_GENES_FILE, symbol_column="gene_symbol") # New gene list

# --- Venn Diagram 1 (Original): YAF Promoter Genes vs. SOX2 Target Genes ---
promoter_overlapping_genes = yaf_promoter_genes.intersection(sox2_genes)
promoter_yaf_only_genes = yaf_promoter_genes.difference(sox2_genes)
promoter_sox2_only_genes = sox2_genes.difference(yaf_promoter_genes)

# --- Venn Diagram 2 (Original): All YAF-Enriched Genes vs. SOX2 Target Genes ---
all_enriched_overlapping_genes = yaf_all_enriched_genes.intersection(sox2_genes)
all_enriched_yaf_only_genes = yaf_all_enriched_genes.difference(sox2_genes)
all_enriched_sox2_only_genes = sox2_genes.difference(yaf_all_enriched_genes)

# Print summary statistics
print(f"\nBasic Statistics:")
print(f"Total All YAF-Enriched genes (H2AK119Ub-YAF2 gene bodies): {len(yaf_all_enriched_genes)}")
print(f"Total YAF Promoter genes (H2AK119Ub-YAF2 promoters): {len(yaf_promoter_genes)}")
print(f"Total SOX2 target genes: {len(sox2_genes)}")
print(f"Total YAF Downregulated genes: {len(yaf_downregulated_genes)}")

print(f"\nOverlap 1 (Original): YAF Promoter Genes vs. SOX2 Target Genes")
print(f"  Number of overlapping genes (YAF Promoter & SOX2): {len(promoter_overlapping_genes)}")
print(f"  Number of YAF Promoter only genes: {len(promoter_yaf_only_genes)}")
print(f"  Number of SOX2 only genes (vs YAF Promoter): {len(promoter_sox2_only_genes)}")

print(f"\nOverlap 2 (Original): All YAF-Enriched Genes vs. SOX2 Target Genes")
print(f"  Number of overlapping genes (All YAF-Enriched & SOX2): {len(all_enriched_overlapping_genes)}")
print(f"  Number of All YAF-Enriched only genes: {len(all_enriched_yaf_only_genes)}")
print(f"  Number of SOX2 only genes (vs All YAF-Enriched): {len(all_enriched_sox2_only_genes)}")

# Create and save Original Venn diagrams
plt.figure(figsize=(16, 8))
plt.subplot(1, 2, 1)
venn2(
    (len(promoter_yaf_only_genes), len(promoter_sox2_only_genes), len(promoter_overlapping_genes)),
    set_labels=('YAF Promoter Genes', 'SOX2 Target Genes'),
    set_colors=('#F9A252', '#64A0C8'), alpha=0.7
)
plt.title('YAF Promoter Genes vs. SOX2 Target Genes', fontsize=14)

plt.subplot(1, 2, 2)
venn2(
    (len(all_enriched_yaf_only_genes), len(all_enriched_sox2_only_genes), len(all_enriched_overlapping_genes)),
    set_labels=('All YAF-Enriched Genes', 'SOX2 Target Genes'),
    set_colors=('#F9A252', '#64A0C8'), alpha=0.7
)
plt.title('All YAF-Enriched Genes vs. SOX2 Target Genes', fontsize=14)

plt.tight_layout(rect=(0, 0, 1, 0.96))
plt.suptitle("YAF and SOX2 Gene Overlap Analysis (2-Way)", fontsize=16, y=0.99)
plt.savefig(os.path.join(OUTPUT_DIR, 'venn_diagrams_2way.png'))
plt.close()

# --- New Venn Diagram 3: All YAF-Enriched vs. SOX2 vs. YAF Downregulated ---
s1_all = yaf_all_enriched_genes
s2_sox2 = sox2_genes
s3_yaf_down = yaf_downregulated_genes

v3_all_s1_only = len(s1_all - (s2_sox2 | s3_yaf_down))
v3_all_s2_only = len(s2_sox2 - (s1_all | s3_yaf_down))
v3_all_s3_only = len(s3_yaf_down - (s1_all | s2_sox2))
v3_all_s1s2 = len((s1_all & s2_sox2) - s3_yaf_down)
v3_all_s1s3 = len((s1_all & s3_yaf_down) - s2_sox2)
v3_all_s2s3 = len((s2_sox2 & s3_yaf_down) - s1_all)
v3_all_s1s2s3 = len(s1_all & s2_sox2 & s3_yaf_down)

# --- New Venn Diagram 4: YAF Promoter vs. SOX2 vs. YAF Downregulated ---
s1_prom = yaf_promoter_genes
# s2_sox2 and s3_yaf_down are the same

v3_prom_s1_only = len(s1_prom - (s2_sox2 | s3_yaf_down))
v3_prom_s2_only = len(s2_sox2 - (s1_prom | s3_yaf_down))
v3_prom_s3_only = len(s3_yaf_down - (s1_prom | s2_sox2))
v3_prom_s1s2 = len((s1_prom & s2_sox2) - s3_yaf_down)
v3_prom_s1s3 = len((s1_prom & s3_yaf_down) - s2_sox2)
v3_prom_s2s3 = len((s2_sox2 & s3_yaf_down) - s1_prom)
v3_prom_s1s2s3 = len(s1_prom & s2_sox2 & s3_yaf_down)

print(f"\nOverlap 3 (New): All YAF-Enriched (A) vs. SOX2 (B) vs. YAF Downregulated (C)")
print(f"  A only (All YAF-Enriched only): {v3_all_s1_only}")
print(f"  B only (SOX2 only): {v3_all_s2_only}")
print(f"  C only (YAF Downregulated only): {v3_all_s3_only}")
print(f"  A & B (not C): {v3_all_s1s2}")
print(f"  A & C (not B): {v3_all_s1s3}")
print(f"  B & C (not A): {v3_all_s2s3}")
print(f"  A & B & C (Intersection): {v3_all_s1s2s3}")

print(f"\nOverlap 4 (New): YAF Promoter (A) vs. SOX2 (B) vs. YAF Downregulated (C)")
print(f"  A only (YAF Promoter only): {v3_prom_s1_only}")
print(f"  B only (SOX2 only): {v3_prom_s2_only}")
print(f"  C only (YAF Downregulated only): {v3_prom_s3_only}")
print(f"  A & B (not C): {v3_prom_s1s2}")
print(f"  A & C (not B): {v3_prom_s1s3}")
print(f"  B & C (not A): {v3_prom_s2s3}")
print(f"  A & B & C (Intersection): {v3_prom_s1s2s3}")

# Create and save New 3-Way Venn diagrams
plt.figure(figsize=(16, 8))
plt.subplot(1, 2, 1)
venn3(
    (v3_all_s1_only, v3_all_s2_only, v3_all_s1s2, v3_all_s3_only, v3_all_s1s3, v3_all_s2s3, v3_all_s1s2s3),
    set_labels=('All YAF-Enriched', 'SOX2 Targets', 'YAF Downregulated'),
    set_colors=('#F9A252', '#64A0C8', '#77DD77')
)
plt.title('All YAF-Enriched vs SOX2 vs YAF Downregulated', fontsize=12)

plt.subplot(1, 2, 2)
venn3(
    (v3_prom_s1_only, v3_prom_s2_only, v3_prom_s1s2, v3_prom_s3_only, v3_prom_s1s3, v3_prom_s2s3, v3_prom_s1s2s3),
    set_labels=('YAF Promoter', 'SOX2 Targets', 'YAF Downregulated'),
    set_colors=('#F9A252', '#64A0C8', '#77DD77')
)
plt.title('YAF Promoter vs SOX2 vs YAF Downregulated', fontsize=12)

plt.tight_layout(rect=(0, 0, 1, 0.96))
plt.suptitle("YAF, SOX2, and YAF Downregulated Gene Overlap Analysis (3-Way)", fontsize=16, y=0.99)
plt.savefig(os.path.join(OUTPUT_DIR, 'venn_diagrams_3way.png'))
plt.close()

# Save gene lists
gene_lists_to_save = {
    'promoter_overlapping_YAF_SOX2.txt': promoter_overlapping_genes,
    'promoter_yaf_only_vs_SOX2.txt': promoter_yaf_only_genes,
    'all_enriched_overlapping_YAF_SOX2.txt': all_enriched_overlapping_genes,
    'all_enriched_yaf_only_vs_SOX2.txt': all_enriched_yaf_only_genes,

    # Lists for "All YAF-Enriched vs SOX2 vs YAF Downregulated"
    'all_enriched_Sox2_YAF_down_INTERSECTION.txt': s1_all & s2_sox2 & s3_yaf_down,
    'all_enriched_ONLY_vs_Sox2_YAF_down.txt': s1_all - (s2_sox2 | s3_yaf_down),
    'Sox2_ONLY_vs_all_enriched_YAF_down.txt': s2_sox2 - (s1_all | s3_yaf_down),
    'YAF_down_ONLY_vs_all_enriched_Sox2.txt': s3_yaf_down - (s1_all | s2_sox2),
    
    # Lists for "YAF Promoter vs SOX2 vs YAF Downregulated"
    'promoter_Sox2_YAF_down_INTERSECTION.txt': s1_prom & s2_sox2 & s3_yaf_down,
    'promoter_ONLY_vs_Sox2_YAF_down.txt': s1_prom - (s2_sox2 | s3_yaf_down),
    'Sox2_ONLY_vs_promoter_YAF_down.txt': s2_sox2 - (s1_prom | s3_yaf_down),
    'YAF_down_ONLY_vs_promoter_Sox2.txt': s3_yaf_down - (s1_prom | s2_sox2),
}

for filename, gene_set in gene_lists_to_save.items():
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        with open(filepath, 'w') as f:
            f.write('\n'.join(sorted(list(gene_set)))) # Ensure it's a list for sorted()
    except Exception as e:
        print(f"Error writing to file {filepath}: {e}")
        sys.exit(1)

def calculate_enrichment_scores():
    """
    Calculate and visualize enrichment scores for regulatory regions.
    Compares H2AK119Ub fold changes between SOX2 targets and non-targets.
    (This function is kept from the original script and might not need changes for the current request)
    """
    try:
        df = pd.read_csv(YAF_FULL_FILE)
        regulatory_regions = ['Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', '5\' UTR', 'Exon']
        first_exon = df['annotation'].str.contains('exon 1 of', na=False)
        regulatory_mask = df['annotation'].isin(regulatory_regions) | first_exon
        
        regulatory_df = df[regulatory_mask].copy()
        regulatory_df['is_sox2_target'] = regulatory_df['SYMBOL'].isin(sox2_genes)
        
        enrichment_stats = regulatory_df.groupby('is_sox2_target')['fold_change'].agg(['mean', 'std', 'count']).round(3)
        enrichment_stats_path = os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_stats.csv')
        enrichment_stats.to_csv(enrichment_stats_path)
        print("\nEnrichment statistics for regulatory regions:")
        print(enrichment_stats)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        sns.boxplot(data=regulatory_df, x='is_sox2_target', y='fold_change', ax=ax1)
        ax1.set_title('With outliers')
        ax1.set_xlabel('Is SOX2 target')
        ax1.set_ylabel('Fold change')
        
        sns.boxplot(data=regulatory_df, x='is_sox2_target', y='fold_change', showfliers=False, ax=ax2)
        ax2.set_title('Without outliers')
        ax2.set_xlabel('Is SOX2 target')
        ax2.set_ylabel('Fold change')
        
        plt.suptitle('H2AK119Ub enrichment in regulatory regions\nSOX2 targets vs non-targets', y=1.05)
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

calculate_enrichment_scores() # Kept from original

print("\nAnalysis complete. Files generated in", OUTPUT_DIR + ":")
print("1. venn_diagrams_2way.png - Visualization of 2-way gene overlaps (YAF vs SOX2)")
print("2. venn_diagrams_3way.png - Visualization of 3-way gene overlaps (YAF vs SOX2 vs YAF Downregulated)")
print("3. Gene lists (*.txt) - Lists of genes in different categories (ready for GO analysis)")
print("4. regulatory_regions_enrichment_stats.csv - Enrichment statistics for regulatory regions")
print("5. regulatory_regions_enrichment_boxplot.png - Visualization of enrichment scores")