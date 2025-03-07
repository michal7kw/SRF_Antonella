#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import glob
import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze overlap between differentially expressed genes and Cut&Tag peaks.')
    parser.add_argument('--deseq2_file', type=str, default='SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv',
                        help='Path to DESeq2 differential expression results')
    parser.add_argument('--counts_file', type=str, default='SRF_RNA/results/deseq2/YAF_vs_GFP/normalized_counts.csv',
                        help='Path to normalized counts file')
    parser.add_argument('--peaks_dir', type=str, 
                        default='/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling',
                        help='Directory containing peak files')
    parser.add_argument('--gtf_file', type=str, default='/beegfs/scratch/ric.broccoli/kubacki.michal/genomes/human/gencode.v38.annotation.gtf',
                        help='Path to GTF annotation file')
    parser.add_argument('--padj_threshold', type=float, default=0.05,
                        help='Adjusted p-value threshold for significant genes')
    parser.add_argument('--log2fc_threshold', type=float, default=1.0,
                        help='Log2 fold change threshold for significant genes')
    parser.add_argument('--output_dir', type=str, default='overlap_results',
                        help='Directory to save results')
    
    return parser.parse_args()

def load_deseq2_results(file_path, padj_threshold, log2fc_threshold):
    """Load DESeq2 results and filter for significant genes."""
    print(f"Loading differential expression data from {file_path}")
    df = pd.read_csv(file_path, index_col=0)
    
    # Filter for significant genes
    up_genes = df[(df['padj'] < padj_threshold) & (df['log2FoldChange'] > log2fc_threshold)]
    down_genes = df[(df['padj'] < padj_threshold) & (df['log2FoldChange'] < -log2fc_threshold)]
    
    print(f"Found {len(up_genes)} upregulated and {len(down_genes)} downregulated genes")
    return df, up_genes, down_genes

def load_normalized_counts(file_path):
    """Load normalized counts data."""
    print(f"Loading normalized counts from {file_path}")
    return pd.read_csv(file_path, index_col=0)

def extract_gene_coordinates(gtf_file, gene_ids):
    """Extract gene coordinates from GTF file for the given gene IDs."""
    print(f"Extracting coordinates for {len(gene_ids)} genes from {gtf_file}")
    
    # Create a temporary file with gene IDs
    with open('temp_gene_ids.txt', 'w') as f:
        for gene_id in gene_ids:
            f.write(f"{gene_id}\n")
    
    # Use grep to extract relevant lines from GTF file (much faster than parsing the whole file)
    os.system(f"grep -f temp_gene_ids.txt {gtf_file} | grep 'gene' > temp_genes.gtf")
    
    # Parse the filtered GTF file
    gene_coords = []
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
                    gene_coords.append({
                        'gene_id': gene_id,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    })
    
    # Clean up temporary files
    os.system('rm temp_gene_ids.txt temp_genes.gtf')
    
    # Convert to DataFrame
    gene_coords_df = pd.DataFrame(gene_coords)
    print(f"Found coordinates for {len(gene_coords_df)} genes")
    return gene_coords_df

def load_peak_files(peaks_dir):
    """Load all peak files from the specified directory."""
    peak_files = glob.glob(os.path.join(peaks_dir, '*_broad_peaks_final.broadPeak'))
    print(f"Found {len(peak_files)} peak files")
    
    peaks_by_sample = {}
    for peak_file in peak_files:
        sample_name = os.path.basename(peak_file).split('_broad_peaks')[0]
        peaks = BedTool(peak_file)
        peaks_by_sample[sample_name] = peaks
        print(f"  {sample_name}: {len(peaks)} peaks")
    
    return peaks_by_sample

def find_overlaps(gene_coords_df, peaks_by_sample):
    """Find overlaps between genes and peaks."""
    print("Finding overlaps between genes and peaks")
    
    # Convert gene coordinates to BedTool
    gene_bed_data = []
    for _, row in gene_coords_df.iterrows():
        gene_bed_data.append([
            row['chrom'], 
            row['start'], 
            row['end'], 
            row['gene_id'], 
            '0', 
            row['strand']
        ])
    
    gene_bed = BedTool(gene_bed_data)
    
    # Find overlaps for each sample
    overlaps_by_sample = {}
    for sample, peaks in peaks_by_sample.items():
        overlaps = gene_bed.intersect(peaks, wa=True, wb=True)
        overlapping_genes = set()
        for overlap in overlaps:
            overlapping_genes.add(overlap[3])  # gene_id is in the 4th field (0-based)
        
        overlaps_by_sample[sample] = overlapping_genes
        print(f"  {sample}: {len(overlapping_genes)} overlapping genes")
    
    return overlaps_by_sample

def create_overlap_summary(up_genes, down_genes, overlaps_by_sample):
    """Create a summary of overlaps between DE genes and peaks."""
    print("Creating overlap summary")
    
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
    
    return summary

def plot_overlap_summary(summary, output_dir):
    """Plot summary of overlaps."""
    print("Plotting overlap summary")
    
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
    plt.ylabel('Percent of DE genes with peak overlap')
    plt.title('Overlap between differentially expressed genes and peaks')
    plt.xticks(x, samples, rotation=45)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, 'overlap_summary.png'), dpi=300)
    plt.close()
    
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
    plt.title('Jaccard similarity between overlapping genes (%)')
    plt.tight_layout()
    
    plt.savefig(os.path.join(output_dir, 'sample_similarity.png'), dpi=300)
    plt.close()

def save_results(summary, output_dir):
    """Save results to output directory."""
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
    
    summary_df.to_csv(os.path.join(output_dir, 'overlap_summary.csv'), index=False)
    
    # Save overlapping gene lists
    for sample, data in summary.items():
        with open(os.path.join(output_dir, f'{sample}_up_overlapping_genes.txt'), 'w') as f:
            for gene in data['up_gene_list']:
                f.write(f"{gene}\n")
        
        with open(os.path.join(output_dir, f'{sample}_down_overlapping_genes.txt'), 'w') as f:
            for gene in data['down_gene_list']:
                f.write(f"{gene}\n")

def main():
    """Main function to run the analysis."""
    args = parse_arguments()
    
    # Load and filter DESeq2 results
    deseq2_df, up_genes, down_genes = load_deseq2_results(
        args.deseq2_file, args.padj_threshold, args.log2fc_threshold
    )
    
    # Load normalized counts
    counts_df = load_normalized_counts(args.counts_file)
    
    # Extract gene coordinates
    all_de_genes = set(up_genes.index).union(set(down_genes.index))
    gene_coords_df = extract_gene_coordinates(args.gtf_file, all_de_genes)
    
    # Load peak files
    peaks_by_sample = load_peak_files(args.peaks_dir)
    
    # Find overlaps
    overlaps_by_sample = find_overlaps(gene_coords_df, peaks_by_sample)
    
    # Create overlap summary
    summary = create_overlap_summary(up_genes, down_genes, overlaps_by_sample)
    
    # Plot results
    plot_overlap_summary(summary, args.output_dir)
    
    # Save results
    save_results(summary, args.output_dir)
    
    print("Analysis complete!")

if __name__ == "__main__":
    main() 