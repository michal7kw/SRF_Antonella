#!/bin/bash
#SBATCH --job-name=7_motif_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_motif_analysis.err"
#SBATCH --output="logs/7_motif_analysis.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub/1_iterative_processing/single

# Create necessary directories
mkdir -p results/motifs
mkdir -p logs/motifs

# Run HOMER motif analysis for each sample
for i in {1,2,3}
do
    echo "Running motif analysis for YAF_${i}..."
    findMotifsGenome.pl \
        analysis/peaks/YAF_${i}_peaks.broadPeak \
        ../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        results/motifs/YAF_${i} \
        -size 200 -mask \
        2> logs/motifs/YAF_${i}.log
done

echo "Motif analysis completed"