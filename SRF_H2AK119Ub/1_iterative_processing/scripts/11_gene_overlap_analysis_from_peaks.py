# This script performs overlap analysis between YAF peaks and SOX2 target genes using pre-annotated peak data.

#### Input files: ####
# - sox2_binding.csv: List of SOX2 target genes
# - peak_annotation.csv:
#     Pre-annotated peak file containing gene associations and genomic context
# ---------------------------------------------------------------------------------------------------------
# seqnames","start","end","width","strand","Conc","Conc_GFP","Conc_YAF","Fold","p.value","FDR","annotation","geneChr","geneStart","geneEnd","geneLength","geneStrand","geneId","transcriptId","distanceToTSS"
# "chr6",113119700,113123406,3707,"*",6.89511090568006,4.93362840925132,7.69687840963268,-2.5535222359897,2.37215678850285e-19,9.81384984971516e-15,"Distal Intergenic",6,112988311,113084454,96144,1,"107986636","ENST00000670062.1",131389
# "chr2",205753075,205755602,2528,"*",7.37727236414482,5.63707505655007,8.14337666216095,-2.35000136379224,6.3651210020658e-19,1.31665710488232e-14,"Promoter (<=1kb)",2,205752468,205763806,11339,1,"8828","ENST00000468256.1",607
# "chr3",181565438,181568638,3201,"*",6.40300689448442,4.01330023651189,7.25834189276477,-2.90581408115723,3.90748610304035e-18,5.3885535856294e-14,"Promoter (1-2kb)",3,181563748,181742884,179137,1,"347689","ENST00000498731.6",1690
# "chr10",76350548,76354205,3658,"*",7.38950645788773,5.80440353281342,8.12650013528447,-2.17986571795513,8.51522561446875e-18,8.704714273874e-14,"Intron (ENST00000611255.5/83938, intron 6 of 6)",10,76324401,76557366,232966,1,"83938","ENST00000598708.1",26147
# "chr3",171190768,171192876,2109,"*",7.14042283422677,5.50698292981968,7.88689145771672,-2.22179451704267,1.05203092430374e-17,8.704714273874e-14,"Exon (ENST00000436636.7/23043, exon 6 of 33)",3,171157551,171225711,68161,2,"23043","ENST00000468757.1",32835
# ---------------------------------------------------------------------------------------------------------

# Output files (in analysis/overlap_analysis_from_peaks/):
# 1. venn_diagrams.png:
#     - Two Venn diagrams showing overlap between YAF and SOX2 genes
#     - One for all genes, one for genes in regulatory regions only
# 2. Gene lists for GO analysis:
#     - all_overlapping_genes.txt: Genes found in both YAF and SOX2 sets
#     - all_yaf_only_genes.txt: Genes unique to YAF set
#     - regulatory_overlapping_genes.txt: Genes in regulatory regions found in both sets
#     - regulatory_yaf_only_genes.txt: Genes in regulatory regions unique to YAF
# 3. regulatory_regions_enrichment_stats.csv:
#     Statistics comparing YAF concentration values between SOX2 targets and non-targets
# 4. regulatory_regions_enrichment_boxplot.png:
#     Visualization of YAF concentration values in regulatory regions

import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
from typing import Set
import mygene
import sys

# Define file paths
OUTPUT_DIR = sys.argv[1]
SOX2_GENES_FILE = sys.argv[2]
ANNOTATION_FILE = sys.argv[3]

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def convert_geneid_to_symbol(gene_ids: list) -> dict:
    """
    Convert gene IDs to gene symbols using mygene
    Returns a dictionary mapping gene IDs to symbols
    """
    mg = mygene.MyGeneInfo()
    try:
        # Query mygene for the gene IDs
        results = mg.querymany(gene_ids, scopes='entrezgene', fields='symbol', species='human')
        
        # Create mapping dictionary
        id_to_symbol = {}
        for result in results:
            if 'symbol' in result and 'query' in result:
                id_to_symbol[str(result['query'])] = result['symbol']
        
        return id_to_symbol
    except Exception as e:
        print(f"Error converting gene IDs to symbols: {str(e)}")
        raise

def read_annotated_peaks() -> pd.DataFrame:
    """
    Read pre-annotated peaks from the existing annotation file
    Returns a DataFrame with peak annotations including gene IDs and symbols
    """
    try:
        # Read the existing peak annotation file
        df = pd.read_csv(ANNOTATION_FILE)
        
        # Convert geneId to string to ensure consistent handling
        df['geneId'] = df['geneId'].astype(str)
        
        # Get unique gene IDs and convert to symbols
        unique_gene_ids = df['geneId'].unique().tolist()
        id_to_symbol_map = convert_geneid_to_symbol(unique_gene_ids)
        
        # Add symbol column
        df['SYMBOL'] = df['geneId'].map(id_to_symbol_map)
        
        # Drop rows where symbol conversion failed
        df = df.dropna(subset=['SYMBOL'])
        
        # Use Conc_YAF as the signal value
        df['signalValue'] = df['Conc_YAF']
        
        # Calculate fold change from the existing Fold column
        df['fold_change'] = df['Fold']
        
        # Add regulatory region flag - include promoter regions up to 3kb
        df['is_regulatory'] = df['annotation'].str.contains(
            'Promoter|5\'|UTR|exon 1|Exon.*1 of', 
            case=False
        )
            
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

if not os.path.exists(ANNOTATION_FILE):
    raise FileNotFoundError(f"Peak annotation file not found: {ANNOTATION_FILE}")

# Read pre-annotated peaks
annotated_peaks_df = read_annotated_peaks()

# Get gene sets (convert to strings when creating sets)
yaf_genes = set(str(gene) for gene in annotated_peaks_df['SYMBOL'].unique())
yaf_regulatory_genes = set(str(gene) for gene in get_regulatory_region_genes(annotated_peaks_df))
sox2_genes = set(str(gene) for gene in read_gene_list(SOX2_GENES_FILE))

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
# Use a built-in style instead of seaborn
plt.style.use('default')  # Changed from 'seaborn' to 'default'
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
    
    # Calculate statistics using Conc_YAF values
    enrichment_stats = regulatory_df.groupby('is_sox2_target')['Conc_YAF'].agg([
        ('mean', 'mean'),
        ('median', 'median'),
        ('std', 'std'),
        ('count', 'count')
    ]).round(3)
    
    enrichment_stats.to_csv(os.path.join(OUTPUT_DIR, 'regulatory_regions_enrichment_stats.csv'))
    print("\nEnrichment statistics for regulatory regions:")
    print(enrichment_stats)
    
    # Use default style instead of seaborn
    plt.style.use('default')
    plt.figure(figsize=(8, 6))
    
    # Create enhanced boxplot with corrected parameters
    sns.boxplot(data=regulatory_df, 
                x='is_sox2_target', 
                y='Conc_YAF',
                color='lightgray')
    
    # Add individual points
    sns.stripplot(data=regulatory_df, 
                 x='is_sox2_target', 
                 y='Conc_YAF',
                 color='darkblue', 
                 alpha=0.3, 
                 size=4, 
                 jitter=0.2)
    
    # Customize appearance
    plt.title('YAF Peak Concentration in Regulatory Regions\nSOX2 Targets vs Non-targets',
              fontsize=12, pad=15)
    plt.xlabel('SOX2 Target Status', fontsize=10)
    plt.ylabel('YAF Concentration', fontsize=10)
    
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