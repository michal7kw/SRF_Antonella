"""
This script performs overlap analysis between YAF peaks (present in at least 2 replicates) 
and SOX2 target genes.

Input files:
- ./analysis/peaks/YAF_{1,2,3}_broad_peaks_final.broadPeak: Raw peak files
- ./sox2_binding.csv: List of SOX2 target genes

Output files (in analysis/overlap_analysis_from_peaks/):
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
from typing import Set

# Define file paths
SOX2_GENES_FILE = "./sox2_binding.csv"
OUTPUT_DIR = "analysis/overlap_analysis_from_peaks"
ANNOTATION_FILE = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/annotation_broad/tables/peak_annotation.csv"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)



def read_annotated_peaks() -> pd.DataFrame:
    """
    Read pre-annotated peaks from the existing annotation file
    Returns a DataFrame with peak annotations including gene IDs and genomic context
    """
    try:
        # Read the existing peak annotation file
        df = pd.read_csv(ANNOTATION_FILE)
        
        # Filter for YAF peaks only
        df = df[df['sample'].str.startswith('YAF_')]
        
        # Extract gene IDs as symbols (since actual gene symbols are not in the file)
        df['SYMBOL'] = df['geneId']
        
        # Calculate fold change from signal value (if not present)
        if 'fold_change' not in df.columns:
            df['fold_change'] = df['signalValue']
        
        # Ensure we have all required columns
        required_columns = ['annotation', 'SYMBOL', 'fold_change']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing required columns in annotation file: {missing_columns}")
        
        # Add regulatory region flag
        df['is_regulatory'] = df['annotation'].str.contains('Promoter|5\'|UTR|exon 1', case=False)
            
        return df
    
    except Exception as e:
        print(f"Error reading peak annotation file: {str(e)}")
        raise

def read_gene_list(file_path: str) -> Set[str]:
    """Read a list of genes from a file, one gene per line"""
    try:
        with open(file_path, 'r') as f:
            genes = set(line.strip() for line in f if line.strip())
        
        if not genes:
            raise ValueError(f"No genes found in {file_path}")
            
        return genes
    except Exception as e:
        print(f"Error reading gene list file: {str(e)}")
        raise

def get_regulatory_region_genes(df: pd.DataFrame) -> Set[str]:
    """Extract genes that are enriched in regulatory regions"""
    return set(df[df['is_regulatory']]['SYMBOL'].unique())

# Validate input files
if not os.path.exists(SOX2_GENES_FILE):
    raise FileNotFoundError(f"SOX2 genes file not found: {SOX2_GENES_FILE}")

# Read pre-annotated peaks
annotated_peaks_df = read_annotated_peaks()

# Get gene sets
yaf_genes = set(annotated_peaks_df['SYMBOL'].unique())
yaf_regulatory_genes = get_regulatory_region_genes(annotated_peaks_df)
sox2_genes = read_gene_list(SOX2_GENES_FILE)

# Find overlapping and unique genes
overlapping_genes = yaf_genes.intersection(sox2_genes)
yaf_only_genes = yaf_genes.difference(sox2_genes)

# Find overlapping and unique genes in regulatory regions
regulatory_overlapping_genes = yaf_regulatory_genes.intersection(sox2_genes)
regulatory_yaf_only_genes = yaf_regulatory_genes.difference(sox2_genes)

# Print summary statistics
print(f"\nBasic Statistics:")
print(f"Total YAF peaks: {len(annotated_peaks_df)}")
print(f"Total YAF genes: {len(yaf_genes)}")
print(f"YAF genes in regulatory regions: {len(yaf_regulatory_genes)}")
print(f"Total SOX2 target genes: {len(sox2_genes)}")
print(f"\nAll genes overlap:")
print(f"Number of overlapping genes: {len(overlapping_genes)}")
print(f"Number of YAF-only genes: {len(yaf_only_genes)}")
print(f"\nRegulatory regions overlap:")
print(f"Number of overlapping regulatory genes: {len(regulatory_overlapping_genes)}")
print(f"Number of YAF-only regulatory genes: {len(regulatory_yaf_only_genes)}")

# Create and save Venn diagrams
# Set style for better visualization
plt.style.use('seaborn')
plt.figure(figsize=(15, 7))

# Create Venn diagrams with improved styling
plt.subplot(1, 2, 1)
venn2([yaf_genes, sox2_genes], 
      set_labels=('YAF genes', 'SOX2 target genes'),
      set_colors=('skyblue', 'lightgreen'),
      alpha=0.7)
plt.title('All Genes Overlap', fontsize=12, pad=15)

plt.subplot(1, 2, 2)
venn2([yaf_regulatory_genes, sox2_genes], 
      set_labels=('YAF regulatory genes', 'SOX2 target genes'),
      set_colors=('skyblue', 'lightgreen'),
      alpha=0.7)
plt.title('Regulatory Regions Overlap', fontsize=12, pad=15)

# Adjust layout and save with high DPI
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'venn_diagrams.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# Save gene lists for downstream GO analysis
gene_lists = {
    'all_overlapping_genes.txt': overlapping_genes,
    'all_yaf_only_genes.txt': yaf_only_genes,
    'regulatory_overlapping_genes.txt': regulatory_overlapping_genes,
    'regulatory_yaf_only_genes.txt': regulatory_yaf_only_genes
}

for filename, gene_set in gene_lists.items():
    with open(os.path.join(OUTPUT_DIR, filename), 'w') as f:
        f.write('\n'.join(sorted(gene_set)))

def calculate_enrichment_scores(df: pd.DataFrame):
    """Calculate and visualize enrichment scores for regulatory regions"""
    # Filter for regulatory regions using the is_regulatory flag
    regulatory_df = df[df['is_regulatory']].copy()
    regulatory_df['is_sox2_target'] = regulatory_df['SYMBOL'].isin(sox2_genes)
    
    # Calculate statistics using signal values
    enrichment_stats = regulatory_df.groupby('is_sox2_target')['signalValue'].agg([
        ('mean', 'mean'),
        ('median', 'median'),
        ('std', 'std'),
        ('count', 'count')
    ]).round(3)
    
    enrichment_stats.to_csv(os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_stats.csv'))
    print("\nEnrichment statistics for regulatory regions:")
    print(enrichment_stats)
    
    # Create visualization
    plt.style.use('seaborn')
    plt.figure(figsize=(8, 6))
    
    # Create enhanced boxplot
    sns.boxplot(data=regulatory_df, x='is_sox2_target', y='signalValue',
                palette=['lightgray', 'lightgreen'])
    
    # Add individual points
    sns.stripplot(data=regulatory_df, x='is_sox2_target', y='signalValue',
                  color='darkblue', alpha=0.3, size=4, jitter=0.2)
    
    # Customize appearance
    plt.title('Peak Signal Values in Regulatory Regions\nSOX2 Targets vs Non-targets',
              fontsize=12, pad=15)
    plt.xlabel('SOX2 Target Status', fontsize=10)
    plt.ylabel('Signal Value', fontsize=10)
    
    # Set x-axis labels
    plt.xticks([0, 1], ['Non-target', 'SOX2 Target'])
    
    # Save plot with high resolution
    plt.savefig(os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_boxplot.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

calculate_enrichment_scores(annotated_peaks_df)

# Print summary of generated files
print("\nAnalysis complete. Files generated in", OUTPUT_DIR + ":")
print("1. venn_diagrams.png - Visualization of gene overlaps")
print("2. *_genes.txt - Lists of genes in different categories (ready for GO analysis)")
print("3. regulatory_regions_enrichment_stats.csv - Enrichment statistics for regulatory regions")
print("4. regulatory_regions_enrichment_boxplot.png - Visualization of enrichment scores") 