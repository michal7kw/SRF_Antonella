#!/bin/bash
# This version is adapted for local execution.

#### Input files: ####
# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_symbols.txt:
#     List of genes enriched in YAF samples (gene symbols)
# ---------------------------------------------------------------------------------------------------------
# MIR8071-2
# MIR8071-1
# TBX4
# TBC1D3P1-DHX40P1
# ATP9B
# MST1L
# ---------------------------------------------------------------------------------------------------------

# - analysis/8_annotation_and_enrichment/gene_lists_broad/YAF_enriched_genes_broad_full.csv:
#     Full data for YAF-enriched genes including annotation and fold change
# ---------------------------------------------------------------------------------------------------------
# "ENTREZID","SYMBOL","distanceToTSS","annotation","fold_change"
# "102466889","MIR8071-2",6171,"Exon (ENST00000497397.1/ENST00000497397.1, exon 2 of 3)",1.76251149051725
# "102465871","MIR8071-1",-31123,"Exon (ENST00000497872.4/ENST00000497872.4, exon 1 of 5)",1.63204716725618
# "9496","TBX4",-7081,"Distal Intergenic",1.6301302738487
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

#### Output files (in ./analysis/11_gene_overlap_analysis/ relative to script execution): ####
# 1. venn_diagrams.png:
#     - Two Venn diagrams showing overlap between YAF and SOX2 genes
#     - One for all genes, one for genes in regulatory regions only
# 2. Gene lists for GO analysis:
#     - all_overlapping_genes.txt: Genes found in both YAF and SOX2 sets
#     - all_yaf_only_genes.txt: Genes unique to YAF set
#     - regulatory_overlapping_genes.txt: Genes in regulatory regions found in both sets
#     - regulatory_yaf_only_genes.txt: Genes in regulatory regions unique to YAF
# 3. regulatory_regions_enrichment_stats.csv:
#     Statistics comparing fold changes between SOX2 targets and non-targets
# 4. regulatory_regions_enrichment_boxplot.png:
#     Visualization of enrichment scores in regulatory regions

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Activate conda environment
# Ensure the correct conda environment (e.g., 'snakemake' or one with Python and required packages like pandas, matplotlib, matplotlib_venn)
# is active before running this script.
# Example: conda activate your_python_env
# source /opt/common/tools/ric.cosr/miniconda3/bin/activate # Removed cluster-specific path
# conda activate snakemake # Assuming environment is already active

# Define working directory
WORKDIR="."
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory $WORKDIR"; exit 1; }

# Define output directory
OUTPUT_DIR="analysis/11_gene_overlap_analysis"

# Create necessary directories
log_message "Creating output directories..."
mkdir -p "${OUTPUT_DIR}" logs

# Run python script for gene overlap analysis
log_message "Running gene overlap analysis..."
python scripts/11_gene_overlap_analysis.py "${OUTPUT_DIR}"

log_message "Gene overlap analysis completed"