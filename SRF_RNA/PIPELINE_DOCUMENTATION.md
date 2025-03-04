# RNA-Seq Analysis Pipeline Documentation

## Overview

RNA-Seq analysis pipeline. The pipeline processes raw RNA-Seq paired-end sequencing data, performs quality control, read alignment, quantification, and differential expression analysis to identify genes that are differentially expressed between experimental conditions.

## Experimental Design

### Samples and Conditions

The experiment consists of 9 samples across 3 experimental conditions:

- **Control (C)**: Samples C1, C2, C3
- **GFP**: Samples GFP1, GFP2, GFP3
- **YAF**: Samples YAF1, YAF2, YAF3

### Comparisons

The following differential expression comparisons were performed:

1. GFP vs Control (C)
2. YAF vs Control (C)
3. YAF vs GFP

## Pipeline Components

The pipeline was implemented using Snakemake.
The workflow consists of the following major steps:

### 1. Quality Control

#### FastQC on Raw Reads

- **Tool**: FastQC v0.11.9
- **Parameters**:
  - Threads: 4
- **Purpose**: Assess the quality of raw sequencing data to identify potential issues such as adapter contamination, low quality bases, or other sequencing artifacts.
- **Output**: HTML reports and ZIP archives containing quality metrics for each sample.

### 2. Read Trimming

#### Trimmomatic

- **Tool**: Trimmomatic v0.39
- **Parameters**:
  - Threads: 16
  - ILLUMINACLIP: TruSeq3-PE.fa:2:30:10 (Adapter sequences removal with seed mismatches:palindrome clip threshold:simple clip threshold)
  - LEADING: 3 (Remove leading low quality bases below quality 3)
  - TRAILING: 3 (Remove trailing low quality bases below quality 3)
  - SLIDINGWINDOW: 4:15 (Scan with a 4-base window, cutting when average quality per base drops below 15)
  - MINLEN: 36 (Drop reads below 36 bases long)
- **Purpose**: Remove adapter sequences and low-quality bases from raw reads to improve alignment accuracy.
- **Output**: Trimmed paired and unpaired FASTQ files for each sample.

#### FastQC on Trimmed Reads

- **Tool**: FastQC v0.11.9
- **Parameters**:
  - Threads: 2
- **Purpose**: Verify the quality improvement after trimming.
- **Output**: HTML reports and ZIP archives containing quality metrics for trimmed reads.

### 3. Genome Indexing

#### STAR Index Generation

- **Tool**: STAR v2.7.10a
- **Parameters**:
  - Threads: 32
  - sjdbOverhang: 100 (Splice junction database overhang)
  - limitGenomeGenerateRAM: 60GB
- **Reference Genome**: Homo sapiens GRCh38 primary assembly (Ensembl)
- **Annotation**: GENCODE v43 basic annotation (without "chr" prefix)
- **Purpose**: Create an index of the reference genome for efficient read alignment.
- **Output**: STAR genome index files.

### 4. Read Alignment

#### STAR Alignment

- **Tool**: STAR v2.7.10a
- **Parameters**:
  - Threads: 32
  - outSAMtype: BAM SortedByCoordinate (Output sorted BAM files)
  - outSAMattributes: Standard (Include standard SAM attributes)
  - limitBAMsortRAM: 45GB (Memory limit for BAM sorting)
- **Purpose**: Align trimmed reads to the reference genome, identifying splice junctions and mapping reads to genomic locations.
- **Output**: Sorted BAM files and alignment statistics for each sample.

### 5. Post-Alignment Processing

#### BAM Indexing

- **Tool**: Samtools v1.15
- **Purpose**: Create indices for BAM files to enable efficient random access.
- **Output**: BAI index files for each BAM file.

#### BigWig Generation

- **Tool**: deepTools bamCoverage v3.5.1
- **Parameters**:
  - Threads: 16
  - binSize: 10 (Resolution of the bigWig file)
  - normalizeUsing: RPKM (Reads Per Kilobase of transcript, per Million mapped reads normalization)
  - effectiveGenomeSize: 2,913,022,398 (Effective size of the human genome)
- **Purpose**: Generate normalized coverage tracks for visualization in genome browsers.
- **Output**: BigWig files for each sample.

### 6. Gene Expression Quantification

#### featureCounts

- **Tool**: featureCounts (Subread package) v2.0.1
- **Parameters**:
  - Threads: 32
  - -p (Count fragments instead of reads for paired-end data)
  - -t exon (Count reads mapping to exons)
  - -g gene_id (Assign counts to genes based on gene_id attribute)
- **Annotation**: GENCODE v43 basic annotation (without "chr" prefix)
- **Purpose**: Count the number of reads/fragments mapping to each gene.
- **Output**: A tab-delimited file containing read counts for all samples.

### 7. Differential Expression Analysis

#### DESeq2

- **Tool**: DESeq2 v1.36.0 (R package)
- **Parameters**:
  - Design: ~condition (Model expression as a function of experimental condition)
  - Reference level: Control condition for each comparison
  - Low count filtering: Genes with fewer than 10 total counts across all samples were removed
- **Purpose**: Identify differentially expressed genes between conditions.
- **Output**:
  - CSV files containing differential expression results for each comparison
  - Normalized count matrices
  - MA plots showing log fold changes vs mean expression
  - Volcano plots highlighting significantly differentially expressed genes
  - R data objects (RDS) for further analysis

### 8. Visualization and Quality Assessment

#### PCA and Heatmaps

- **Tool**: Custom R scripts using DESeq2, ggplot2, and pheatmap
- **Visualizations**:
  - PCA plot of all samples based on variance-stabilized transformed data
  - Sample correlation heatmap
  - Heatmap of top differentially expressed genes across all comparisons
- **Purpose**: Assess sample clustering, identify potential outliers, and visualize expression patterns of top differentially expressed genes.
- **Output**: PDF files containing the visualizations.

#### MultiQC

- **Tool**: MultiQC v1.12
- **Purpose**: Aggregate quality control metrics from various steps of the pipeline into a single report.
- **Input**: FastQC reports, STAR alignment logs, featureCounts summary
- **Output**: HTML report containing interactive visualizations of quality metrics.

## Computational Requirements

The pipeline was executed on a high-performance computing cluster with the following resource allocations:

- **STAR Indexing**: 32 CPUs, 64GB RAM, 12 hours
- **STAR Alignment**: 32 CPUs, 48GB RAM, 8 hours per sample
- **Trimmomatic**: 16 CPUs, 16GB RAM, 4 hours per sample
- **featureCounts**: 32 CPUs, 32GB RAM, 6 hours
- **DESeq2 Analysis**: 1 CPU, 32GB RAM, 2 hours per comparison
- **BigWig Generation**: 16 CPUs, 32GB RAM, 4 hours per sample

## File Organization

The pipeline organizes output files in a structured directory hierarchy:

```
SRF_RNA/
├── DATA/                      # Raw input data
│   ├── C1/
│   ├── C2/
│   └── ...
├── results/                   # Analysis results
│   ├── fastqc/                # FastQC reports
│   ├── trimmed/               # Trimmed reads
│   ├── star/                  # STAR alignment results
│   │   ├── C1/
│   │   ├── C2/
│   │   └── ...
│   ├── counts/                # featureCounts results
│   ├── bigwig/                # Coverage tracks
│   ├── deseq2/                # Differential expression results
│   │   ├── GFP_vs_C/
│   │   ├── YAF_vs_C/
│   │   ├── YAF_vs_GFP/
│   │   └── ...
│   └── multiqc/               # MultiQC report
├── logs/                      # Log files
│   ├── fastqc/
│   ├── trimmomatic/
│   └── ...
└── scripts/                   # Analysis scripts
```

## Reference Data

The pipeline uses the following reference data:

- **Genome**: Homo sapiens GRCh38 primary assembly
  - Path: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- **Annotation**: GENCODE v43 basic annotation (without "chr" prefix)
  - Path: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.nochr.gtf`
- **Chromosome Sizes**: Derived from the reference genome
  - Path: `/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes`

## Differential Expression Analysis Details

### Statistical Method

The differential expression analysis was performed using DESeq2, which implements a method based on the negative binomial distribution. The analysis workflow includes:

1. **Count Normalization**: DESeq2 estimates size factors to account for differences in sequencing depth between samples.
2. **Dispersion Estimation**: Gene-wise dispersion estimates are calculated and then shrunk towards a fitted dispersion-mean relation using an empirical Bayes approach.
3. **Model Fitting**: A generalized linear model is fitted for each gene using the negative binomial distribution with the experimental design.
4. **Hypothesis Testing**: Wald tests are performed to identify genes with significant differences between conditions.
5. **Multiple Testing Correction**: P-values are adjusted for multiple testing using the Benjamini-Hochberg procedure to control the false discovery rate.

### Filtering Criteria

For the identification of differentially expressed genes, the following criteria were applied:

- **Adjusted p-value (padj)**: < 0.05
- **Log2 Fold Change**: Absolute value > 1 (corresponding to a 2-fold change)
- **Base Mean**: Genes with very low counts across all samples (sum < 10) were filtered out before testing.

### Visualization

For each comparison, the following visualizations were generated:

- **MA Plots**: Display the log2 fold changes against the mean of normalized counts, highlighting genes that pass the significance threshold.
- **Volcano Plots**: Show the relationship between statistical significance (-log10 adjusted p-value) and magnitude of change (log2 fold change), with significant genes highlighted.
- **PCA Plot**: Visualizes sample clustering based on the overall gene expression profiles.
- **Heatmaps**: Show expression patterns of the top 50 differentially expressed genes from each comparison.

## Running the Pipeline

The pipeline is executed using a SLURM job scheduler with the following command:

```bash
sbatch run_snakefile.sh
```

The script performs the following steps:

1. Creates necessary output directories
2. Checks for the existence of the STAR index and input files
3. Activates the required conda environment
4. Unlocks the Snakemake working directory if necessary
5. Runs Snakemake with appropriate parameters for SLURM job submission

## Software Dependencies

The pipeline requires the following software:

- **Snakemake**: Workflow management system
- **FastQC**: Quality control tool for high throughput sequence data
- **Trimmomatic**: Flexible read trimming tool for Illumina NGS data
- **STAR**: Spliced Transcripts Alignment to a Reference
- **Samtools**: Suite of programs for interacting with high-throughput sequencing data
- **featureCounts**: Read count program
- **deepTools**: Suite of tools for the visualization, quality control, and analysis of high-throughput sequencing data
- **MultiQC**: Aggregate results from bioinformatics analyses into a single report
- **R packages**:
  - DESeq2: Differential gene expression analysis
  - tidyverse: Collection of R packages for data science
  - ggplot2: Data visualization package
  - ggrepel: Text and label geoms for 'ggplot2'
  - pheatmap: Pretty heatmaps
  - RColorBrewer: ColorBrewer palettes
  - EnhancedVolcano: Volcano plots with improved features