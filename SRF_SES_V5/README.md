# CUT&TAG Analysis Pipeline for SES-V5 Samples

This Snakemake pipeline is designed to analyze CUT&TAG data, specifically for SES-V5 samples.

## Important Note on Data Availability

**Data Limitation**: The original dataset included two SES-V5 samples, but due to data corruption:
- SES-V5ChIP-Seq1_S5 is excluded from analysis (corrupted R1 fastq file)
- Only SES-V5ChIP-Seq2_S6 is analyzed as a representative sample
- InputSES-V5ChIP-Seq_S2 is used as the input control

This limitation should be considered when interpreting the results, as we lack biological replication.

## Pipeline Steps

The pipeline includes the following steps:
1. Quality control (FastQC and MultiQC)
2. Alignment to reference genome (Bowtie2)
3. Peak calling (MACS2)
4. BigWig file generation for visualization

## Output Files

The pipeline generates the following output directories:
- `qc/`: FastQC and MultiQC reports
- `alignment/`: Aligned BAM files
- `peaks/`: MACS2 peak calls
- `bigwig/`: BigWig files for visualization