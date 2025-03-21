# YAF Enrichment Analysis

This repository contains a streamlined analysis pipeline to examine the enrichment of YAF vs GFP Cut&Tag data at promoters of differentially expressed genes.

## Overview

The analysis aims to determine if there's a relationship between YAF binding and gene expression changes. It compares binding patterns between upregulated and downregulated genes to understand YAF's potential role as a transcriptional activator or repressor.

## Requirements

The analysis requires the following R packages:
- DESeq2
- rtracklayer
- GenomicRanges
- ggplot2
- dplyr
- tidyr
- pheatmap
- RColorBrewer
- gridExtra
- ggpubr
- scales
- viridis

The script will attempt to install any missing packages.

## Input Data

The analysis uses the following input data:
1. DESeq2 results from RNA-seq analysis: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_RNA/results/deseq2/YAF_vs_GFP/differential_expression.csv`
2. GTF annotation file: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.gtf`
3. BigWig files for YAF and GFP samples: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing/analysis/6_bigwig/`

## How to Run

Simply execute the R script:

```bash
chmod +x yaf_enrichment_analysis.R
./yaf_enrichment_analysis.R
```

Or run it with Rscript:

```bash
Rscript yaf_enrichment_analysis.R
```

## Analysis Steps

1. Load DESeq2 results and filter for significant DEGs
2. Load GTF annotation and extract gene information
3. Define promoter regions (2kb upstream and 500bp downstream of TSS)
4. Map DEGs to promoter regions
5. Extract signal from bigWig files for each promoter region
6. Calculate average signal and log2FC for YAF vs GFP
7. Create visualizations
8. Perform statistical tests
9. Generate summary report

## Output

The analysis creates a directory called `YAF_enrichment_results` with the following outputs:

### Data Files
- `YAF_GFP_binding_at_promoters.csv`: Table with binding data for all genes
- `enrichment_summary_stats.csv`: Summary statistics for each regulation group

### Visualizations (in the `plots` subdirectory)
1. Boxplot of YAF/GFP enrichment by regulation status
2. Scatterplot of YAF vs GFP signal
3. Bar plot of average YAF and GFP signal with error bars
4. Violin plot of YAF/GFP ratio distribution
5. Correlation plot between binding and expression log2FC
6. Heatmap of top genes by YAF binding
7. Density plot of YAF and GFP signal distribution
8. Combined figure with multiple plots

### Reports
- `analysis_summary.txt`: Text summary of the analysis results, including statistical tests and biological interpretation

## Interpretation

The analysis helps determine:
1. Whether YAF binding is enriched at promoters of upregulated or downregulated genes
2. If there's a correlation between YAF binding and gene expression changes
3. The potential role of YAF as a transcriptional activator or repressor

## Troubleshooting

If you encounter issues with missing bigWig files, check that the file paths are correct and that the files have the expected names (YAF_1.bw, YAF_2.bw, YAF_3.bw, GFP_1.bw, GFP_2.bw, GFP_3.bw).

If you get errors related to package installation, you may need to install them manually:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "rtracklayer", "GenomicRanges"))
install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer", "gridExtra", "ggpubr", "scales", "viridis"))
```

# Gene Ontology Analysis for YAF vs GFP

This folder contains scripts to perform Gene Ontology (GO) analysis on differentially expressed genes from the YAF vs GFP comparison.

## Files

- `go_analysis.py`: Main Python script that performs GO analysis on upregulated and downregulated genes
- `requirements.txt`: List of Python package dependencies
- `run_go_analysis.sh`: Shell script to set up the environment and run the analysis

## Workflow

The analysis performs the following steps:

1. Loads differential expression results from the DESeq2
2. Identifies significantly up and down-regulated genes (padj < 0.05)
3. Performs GO term enrichment analysis using the EnrichR API (via gseapy)
4. Creates visualizations of top enriched GO terms
5. Performs Gene Set Enrichment Analysis (GSEA) on the ranked gene list

## Categories of GO Analysis

For both upregulated and downregulated genes, the script performs enrichment for:

- Biological Process (BP)
- Cellular Component (CC)
- Molecular Function (MF)

## Running the Analysis

1. Make sure the shell script is executable:
   ```
   chmod +x run_go_analysis.sh
   ```

2. Run the script:
   ```
   ./run_go_analysis.sh
   ```

3. Results will be saved in the `go_analysis_results` directory:
   - Enrichment results in TSV format
   - Visualization plots in PNG format
   - Lists of up and downregulated genes 