# Gene Promoter-Peak Overlap Analysis

This script analyzes the overlap between differentially expressed genes' promoters from RNA-seq data and peaks from Cut&Tag experiments.

## Requirements

- Python 3.6+
- pandas
- numpy
- matplotlib
- seaborn
- pybedtools

You can install the required packages with:

```bash
pip install pandas numpy matplotlib seaborn pybedtools
```

## Usage

```bash
python overlap_analysis.py [options]
```

### Options

- `--deseq2_file`: Path to DESeq2 differential expression results (default: 'SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv')
- `--counts_file`: Path to normalized counts file (default: 'SRF_RNA/results/deseq2/YAF_vs_GFP/normalized_counts.csv')
- `--peaks_dir`: Directory containing peak files (default: '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/5_peak_calling')
- `--gtf_file`: Path to GTF annotation file (default: '/beegfs/scratch/ric.broccoli/kubacki.michal/genomes/human/gencode.v38.annotation.gtf')
- `--padj_threshold`: Adjusted p-value threshold for significant genes (default: 0.05)
- `--log2fc_threshold`: Log2 fold change threshold for significant genes (default: 1.0)
- `--output_dir`: Directory to save results (default: 'overlap_results')
- `--promoter_upstream`: Distance upstream of TSS to define promoter region (default: 2000)
- `--promoter_downstream`: Distance downstream of TSS to define promoter region (default: 200)

## Example

```bash
python overlap_analysis.py --padj_threshold 0.01 --log2fc_threshold 1.5 --promoter_upstream 3000 --output_dir my_results
```

## Output

The script generates the following outputs in the specified output directory:

1. `promoter_overlap_summary.csv`: A summary of the overlap between differentially expressed genes' promoters and peaks for each sample
2. `promoter_overlap_summary.png`: A bar plot showing the percentage of differentially expressed genes with promoter-peak overlaps
3. `promoter_sample_similarity.png`: A heatmap showing the similarity between samples based on overlapping genes
4. `{sample}_up_promoter_overlapping_genes.txt`: List of upregulated genes with promoter-peak overlaps for each sample
5. `{sample}_down_promoter_overlapping_genes.txt`: List of downregulated genes with promoter-peak overlaps for each sample

## How it works

1. Loads differentially expressed genes from DESeq2 results
2. Extracts promoter coordinates for these genes from a GTF file
   - For genes on the + strand: TSS - upstream to TSS + downstream
   - For genes on the - strand: TES - downstream to TES + upstream
3. Loads peak files for each sample
4. Finds overlaps between promoters and peaks using pybedtools
5. Summarizes and visualizes the results

## Notes

- The script expects broadPeak files with the naming convention `{sample}_broad_peaks_final.broadPeak`
- Gene IDs in the DESeq2 results should match those in the GTF file
- By default, promoters are defined as 2000bp upstream to 200bp downstream of the transcription start site (TSS) 