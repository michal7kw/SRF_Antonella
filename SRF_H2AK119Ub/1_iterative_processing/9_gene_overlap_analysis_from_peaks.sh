#!/bin/bash
#SBATCH --job-name=9_gene_overlap_analysis_from_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/9_gene_overlap_analysis_from_peaks.err"
#SBATCH --output="logs/9_gene_overlap_analysis_from_peaks.out"

# This script performs overlap analysis between YAF peaks and SOX2 target genes using pre-annotated peak data.

#### Input files: ####
# - ./sox2_binding.csv: List of SOX2 target genes
# - ./analysis/annotation_broad/tables/peak_annotation.csv:
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

set -e
set -u
set -o pipefail

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/SRF_H2AK119Ub/1_iterative_processing"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

# Create necessary directories
log_message "Creating output directories..."
mkdir -p analysis/overlap_analysis_from_peaks logs

# Run python script for gene overlap analysis
log_message "Running gene overlap analysis..."
python scripts/9_gene_overlap_analysis_from_peaks.py

log_message "Gene overlap analysis completed" 