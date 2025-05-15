# This script performs overlap analysis between YAF-enriched genes and SOX2 target genes.

#### Input files: ####
# - YAF_enriched_genes_broad_promoters.csv:
#     Data for YAF-enriched gene promoters including annotation, fold change, and statistical significance
# ---------------------------------------------------------------------------------------------------------
# "ENTREZID","SYMBOL","distanceToTSS","annotation","fold_change","p_value","FDR"
# "407047","MIR9-2",0,"Promoter (<=1kb)",2.88394266788565,4.62621671744608e-20,4.13738219019504e-15
# "57142","RTN4",2801,"Promoter (2-3kb)",2.85787059102815,1.23987698289329e-11,5.14809227420024e-09
# "347689","SOX2-OT",1690,"Promoter (1-2kb)",2.81379453300378,9.30098704665971e-18,1.44543982267268e-13
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
#     - Venn diagram showing overlap between YAF promoter binding sites and SOX2 target genes
# 2. Gene lists for GO analysis:
#     - overlapping_genes.txt: Genes found in both YAF and SOX2 sets
#     - yaf_only_genes.txt: Genes unique to YAF set
# 3. enrichment_stats.csv:
#     Statistics comparing fold changes between SOX2 targets and non-targets
# 4. enrichment_boxplot.png:
#     Visualization of enrichment scores

import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')

# Define file paths and create output directory
OUTPUT_DIR = sys.argv[1]

# BASE_PATH = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5"
# BASE_PATH = "D:/Github/SRF_H2AK119Ub_cross_V5"
BASE_PATH = "/mnt/d/Github/SRF_H2AK119Ub_cross_V5"

YAF_PROMOTERS_FILE = os.path.join(BASE_PATH, "SRF_H2AK119Ub/1_iterative_processing/analysis/8_annotation_and_enrichment_strict/gene_lists_broad/YAF_enriched_genes_broad_promoters.csv")
SOX2_GENES_FILE = os.path.join(BASE_PATH, "COMMON_DATA/sox2_binding.csv")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def read_gene_list(file_path):
    """Read a list of genes from a file, one gene per line"""
    try:
        with open(file_path, 'r', encoding='latin-1') as f:
            return set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def get_significant_promoter_genes(file_path, fdr_threshold=0.05):
    """
    Extract genes that have statistically significant enrichment in promoter regions
    based on the specified FDR threshold.
    
    Returns a set of gene symbols
    """
    try:
        df = pd.read_csv(file_path)
        # Filter for statistically significant sites
        significant_sites = df[df['FDR'] <= fdr_threshold]
        
        # Extract gene symbols
        return set(significant_sites['SYMBOL'].unique())
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

# Read gene lists and find overlaps
sox2_genes = read_gene_list(SOX2_GENES_FILE)
yaf_promoter_genes = get_significant_promoter_genes(YAF_PROMOTERS_FILE)

# Find overlapping and unique genes for promoter regions
overlapping_genes = yaf_promoter_genes.intersection(sox2_genes)
yaf_only_genes = yaf_promoter_genes.difference(sox2_genes)

# Print summary statistics
print(f"\nBasic Statistics:")
print(f"YAF genes with significant promoter binding: {len(yaf_promoter_genes)}")
print(f"Total SOX2 target genes: {len(sox2_genes)}")
print(f"\nOverlap analysis:")
print(f"Number of overlapping genes: {len(overlapping_genes)}")
print(f"Number of YAF-only genes: {len(yaf_only_genes)}")

# Create and save Venn diagram
plt.figure(figsize=(8, 6))
venn2([yaf_promoter_genes, sox2_genes], ('YAF promoter genes', 'SOX2 target genes'))
plt.title('Overlap between YAF promoter binding sites and SOX2 targets')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'venn_diagram.png'))
plt.close()

# Save gene lists for downstream GO analysis
gene_lists = {
    'overlapping_genes.txt': overlapping_genes,
    'yaf_only_genes.txt': yaf_only_genes
}

for filename, gene_set in gene_lists.items():
    filepath = os.path.join(OUTPUT_DIR, filename)
    try:
        with open(filepath, 'w') as f:
            f.write('\n'.join(sorted(gene_set)))
    except Exception as e:
        print(f"Error writing to file {filepath}: {e}")
        sys.exit(1)

def calculate_enrichment_scores(fdr_threshold=0.05):
    """
    Calculate and visualize enrichment scores for promoter regions.
    Compares H2AK119Ub fold changes between SOX2 targets and non-targets.
    
    Outputs:
    - CSV file with mean, std, and count statistics
    - Two boxplots visualization of fold changes (with and without outliers)
    """
    try:
        df = pd.read_csv(YAF_PROMOTERS_FILE)
        
        # Filter for statistically significant sites
        significant_sites = df[df['FDR'] <= fdr_threshold].copy()
        
        # Add SOX2 target information
        significant_sites['is_sox2_target'] = significant_sites['SYMBOL'].isin(sox2_genes)
        
        # Calculate statistics
        enrichment_stats = significant_sites.groupby('is_sox2_target')['fold_change'].agg(['mean', 'std', 'count']).round(3)
        enrichment_stats_path = os.path.join(OUTPUT_DIR, 'enrichment_stats.csv')
        enrichment_stats.to_csv(enrichment_stats_path)
        print("\nEnrichment statistics for promoter binding sites:")
        print(enrichment_stats)
        
        # Create visualization with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # First plot - with outliers
        sns.boxplot(data=significant_sites, x='is_sox2_target', y='fold_change', ax=ax1)
        ax1.set_title('With outliers')
        ax1.set_xlabel('Is SOX2 target')
        ax1.set_ylabel('Fold change')
        
        # Second plot - without outliers
        sns.boxplot(data=significant_sites, x='is_sox2_target', y='fold_change', 
                    showfliers=False, ax=ax2)
        ax2.set_title('Without outliers')
        ax2.set_xlabel('Is SOX2 target')
        ax2.set_ylabel('Fold change')
        
        # Add overall title
        plt.suptitle('H2AK119Ub enrichment in promoter regions\nSOX2 targets vs non-targets', 
                     y=1.05)
        
        plt.tight_layout()
        enrichment_boxplot_path = os.path.join(OUTPUT_DIR, 'enrichment_boxplot.png')
        plt.savefig(enrichment_boxplot_path, bbox_inches='tight')
        plt.close()
    
    except FileNotFoundError:
        print(f"Error: The file {YAF_PROMOTERS_FILE} was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during enrichment score calculation: {e}")
        sys.exit(1)

calculate_enrichment_scores()

# Print summary of generated files
print("\nAnalysis complete. Files generated in", OUTPUT_DIR + ":")
print("1. venn_diagram.png - Visualization of gene overlaps")
print("2. overlapping_genes.txt and yaf_only_genes.txt - Lists of genes (ready for GO analysis)")
print("3. enrichment_stats.csv - Enrichment statistics for promoter regions")
print("4. enrichment_boxplot.png - Visualization of enrichment scores")
