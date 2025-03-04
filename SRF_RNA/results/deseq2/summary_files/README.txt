DESeq2 Analysis Summary Files
===========================

This directory contains summary files from the DESeq2 differential expression analysis.

File Descriptions:
----------------

master_summary.csv - Overview of all comparisons with key statistics
  - comparison: Name of the comparison
  - total_genes: Total number of genes analyzed
  - up_regulated: Number of significantly up-regulated genes
  - down_regulated: Number of significantly down-regulated genes
  - significant_genes: Total number of significantly differentially expressed genes
  - max_log2fc: Maximum log2 fold change among significant genes
  - min_log2fc: Minimum log2 fold change among significant genes
  - top_up_gene: Gene symbol with highest positive fold change
  - top_down_gene: Gene symbol with highest negative fold change

For each comparison, the following files are generated:

1. [comparison]_detailed_results.csv - Complete DESeq2 results for all genes
   - gene_id: Original gene identifier
   - gene_symbol: Corresponding gene symbol
   - baseMean: Average normalized count across all samples
   - log2FoldChange: Log2 fold change
   - lfcSE: Standard error of log2 fold change
   - stat: Wald statistic
   - pvalue: Raw p-value
   - padj: Adjusted p-value (FDR)

2. [comparison]_significant_genes.csv - Genes with adjusted p-value < cutoff

3. [comparison]_up_regulated.csv - Significantly up-regulated genes

4. [comparison]_down_regulated.csv - Significantly down-regulated genes

Analysis performed with cutoffs: adjusted p-value < 0.05 and absolute log2 fold change > 1
Date: 2025-03-04
