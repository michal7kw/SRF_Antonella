# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import glob
import argparse


# %%
file_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv"
padj_threshold = 0.01
log2fc_threshold = 1.5
output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/overlap_Ub"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(file_path, index_col=0)

# Filter for significant genes
up_genes = df[(df['padj'] < padj_threshold) & (df['log2FoldChange'] > log2fc_threshold)]
down_genes = df[(df['padj'] < padj_threshold) & (df['log2FoldChange'] < -log2fc_threshold)]

print(f"Found {len(up_genes)} upregulated and {len(down_genes)} downregulated genes")

# %%
up_genes.head()

# %%
down_genes.head()

# %%
counts_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/normalized_counts.csv"
counts_df = pd.read_csv(counts_file, index_col=0)

# %%
counts_df.head()

# %%
all_de_genes = set(up_genes.index).union(set(down_genes.index))

# %%
gtf_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf"
promoter_upstream=2000
promoter_downstream=200
gene_ids = all_de_genes

# Create a temporary file with gene IDs
with open('temp_gene_ids.txt', 'w') as f:
    for gene_id in gene_ids:
        f.write(f"{gene_id}\n")

# Use grep to extract relevant lines from GTF file (much faster than parsing the whole file)
os.system(f"grep -f temp_gene_ids.txt {gtf_file} | grep 'gene' > temp_genes.gtf")

# Parse the filtered GTF file to get promoter regions
promoter_coords = []
with open('temp_genes.gtf', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if fields[2] == 'gene':
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            
            # Extract gene_id from attributes
            attributes = fields[8]
            gene_id = None
            for attr in attributes.split(';'):
                if 'gene_id' in attr:
                    gene_id = attr.split('"')[1]
                    break
            
            if gene_id and gene_id in gene_ids:
                # Define promoter region based on strand
                if strand == '+':
                    promoter_start = max(1, start - promoter_upstream)
                    promoter_end = start + promoter_downstream
                else:  # strand == '-'
                    promoter_start = max(1, end - promoter_downstream)
                    promoter_end = end + promoter_upstream
                
                promoter_coords.append({
                    'gene_id': gene_id,
                    'chrom': chrom,
                    'start': promoter_start,
                    'end': promoter_end,
                    'strand': strand
                })

# Clean up temporary files
os.system('rm temp_gene_ids.txt temp_genes.gtf')

# Convert to DataFrame
promoter_coords_df = pd.DataFrame(promoter_coords)
print(f"Found promoter coordinates for {len(promoter_coords_df)} genes")

# %%
promoter_coords_df.head()

# %%
# peaks_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/peaks_Ub"
peaks_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling"
peak_files = glob.glob(os.path.join(peaks_dir, '*_broad_peaks_final.broadPeak'))
print(f"Found {len(peak_files)} peak files")

peaks_by_sample = {}
for peak_file in peak_files:
    sample_name = os.path.basename(peak_file).split('_broad_peaks')[0]
    try:
        # Load the peak file into a pandas DataFrame
        peak_df = pd.read_csv(peak_file, sep='\t', header=None, low_memory=False)

        # Convert the start and end columns to integers
        peak_df[1] = peak_df[1].astype(int)
        peak_df[2] = peak_df[2].astype(int)

        # Create a BedTool from the DataFrame
        peaks = BedTool.from_dataframe(peak_df)
        
        peaks_by_sample[sample_name] = peaks
        print(f"  {sample_name}: {len(peaks)} peaks")
    except Exception as e:
        print(f"Error processing {peak_file}: {e}")
        print(f"Error details: {type(e).__name__}, {e}")
        print("Skipping this file.")

# %%
# First, let's standardize the chromosome names in the promoter data
promoter_bed_data = []
for _, row in promoter_coords_df.iterrows():
    chrom = str(row['chrom'])
    # Remove 'chr' prefix if present to standardize naming
    chrom = chrom.replace('chr', '')
    
    promoter_bed_data.append([
        chrom,               # Standardized chromosome name
        int(row['start']),   
        int(row['end']),    
        str(row['gene_id']), 
        '0',
        str(row['strand'])   
    ])

promoter_bed = BedTool(promoter_bed_data)

# Find overlaps for each sample
overlaps_by_sample = {}
for sample, peaks in peaks_by_sample.items():
    try:
        # Convert peaks to a new BedTool with standardized chromosome names
        peaks_data = []
        for peak in peaks:
            chrom = str(peak[0]).replace('chr', '')  # Standardize chromosome names in peaks
            peaks_data.append([chrom] + list(peak[1:]))
        
        standardized_peaks = BedTool(peaks_data)
        
        overlaps = promoter_bed.intersect(standardized_peaks, wa=True, wb=True)
        overlapping_genes = set()
        for overlap in overlaps:
            overlapping_genes.add(overlap[3])  # gene_id is in the 4th field (0-based)
        
        overlaps_by_sample[sample] = overlapping_genes
        print(f"  {sample}: {len(overlapping_genes)} genes with promoter overlaps")
    except Exception as e:
        print(f"Error during intersection for sample {sample}: {e}")

# %%
up_gene_ids = set(up_genes.index)
down_gene_ids = set(down_genes.index)

summary = {}
for sample, overlapping_genes in overlaps_by_sample.items():
    up_overlaps = up_gene_ids.intersection(overlapping_genes)
    down_overlaps = down_gene_ids.intersection(overlapping_genes)
    
    summary[sample] = {
        'up_overlaps': len(up_overlaps),
        'up_genes': len(up_gene_ids),
        'up_percent': len(up_overlaps) / len(up_gene_ids) * 100 if len(up_gene_ids) > 0 else 0,
        'down_overlaps': len(down_overlaps),
        'down_genes': len(down_gene_ids),
        'down_percent': len(down_overlaps) / len(down_gene_ids) * 100 if len(down_gene_ids) > 0 else 0,
        'up_gene_list': list(up_overlaps),
        'down_gene_list': list(down_overlaps)
    }

# %%
# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Prepare data for plotting
samples = list(summary.keys())
up_percents = [summary[s]['up_percent'] for s in samples]
down_percents = [summary[s]['down_percent'] for s in samples]

# Create bar plot
plt.figure(figsize=(12, 6))
x = np.arange(len(samples))
width = 0.35

plt.bar(x - width/2, up_percents, width, label='Upregulated genes')
plt.bar(x + width/2, down_percents, width, label='Downregulated genes')

plt.xlabel('Sample')
plt.ylabel('Percent of DE genes with promoter peak overlap')
plt.title('Overlap between differentially expressed genes promoters and peaks')
plt.xticks(x, samples, rotation=45)
plt.legend()
plt.tight_layout()

plt.show()

# Create heatmap for sample comparison
overlap_matrix = np.zeros((len(samples), len(samples)))
for i, sample1 in enumerate(samples):
    for j, sample2 in enumerate(samples):
        if i == j:
            overlap_matrix[i, j] = 100
        else:
            genes1 = set(summary[sample1]['up_gene_list'] + summary[sample1]['down_gene_list'])
            genes2 = set(summary[sample2]['up_gene_list'] + summary[sample2]['down_gene_list'])
            
            if len(genes1) > 0 and len(genes2) > 0:
                overlap = len(genes1.intersection(genes2))
                union = len(genes1.union(genes2))
                overlap_matrix[i, j] = overlap / union * 100

plt.figure(figsize=(10, 8))
sns.heatmap(overlap_matrix, annot=True, fmt='.1f', cmap='viridis',
            xticklabels=samples, yticklabels=samples)
plt.title('Jaccard similarity between samples based on promoter overlaps (%)')
plt.tight_layout()

plt.show()

# %%
print(f"Saving results to {output_dir}")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Save summary as CSV
summary_df = pd.DataFrame({
    'Sample': [],
    'Upregulated_Overlaps': [],
    'Upregulated_Total': [],
    'Upregulated_Percent': [],
    'Downregulated_Overlaps': [],
    'Downregulated_Total': [],
    'Downregulated_Percent': []
})

for sample, data in summary.items():
    summary_df = pd.concat([summary_df, pd.DataFrame({
        'Sample': [sample],
        'Upregulated_Overlaps': [data['up_overlaps']],
        'Upregulated_Total': [data['up_genes']],
        'Upregulated_Percent': [data['up_percent']],
        'Downregulated_Overlaps': [data['down_overlaps']],
        'Downregulated_Total': [data['down_genes']],
        'Downregulated_Percent': [data['down_percent']]
    })], ignore_index=True)

summary_df.to_csv(os.path.join(output_dir, 'promoter_overlap_summary.csv'), index=False)

# Save overlapping gene lists
for sample, data in summary.items():
    with open(os.path.join(output_dir, f'{sample}_up_promoter_overlapping_genes.txt'), 'w') as f:
        for gene in data['up_gene_list']:
            f.write(f"{gene}\n")
    
    with open(os.path.join(output_dir, f'{sample}_down_promoter_overlapping_genes.txt'), 'w') as f:
        for gene in data['down_gene_list']:
            f.write(f"{gene}\n")

# %%



