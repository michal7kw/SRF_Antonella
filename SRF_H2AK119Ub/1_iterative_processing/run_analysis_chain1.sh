#!/bin/bash
#SBATCH --job-name=run_all_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/run_all_analysis.err"
#SBATCH --output="logs/run_all_analysis.out"

# Set error handling
set -e  # Exit on error
set -u  # Exit on undefined variable

echo "Submitting analysis pipeline at $(date)"

# Submit first job
job1=$(sbatch 5_differential_binding.sh | awk '{print $4}')
echo "Submitted 5_differential_binding.sh (JobID: $job1)"

# Submit subsequent jobs with dependencies
job2=$(sbatch --dependency=afterok:$job1 6_annotation_and_enrichment.sh | awk '{print $4}')
echo "Submitted 6_annotation_and_enrichment.sh (JobID: $job2)"

job3=$(sbatch --dependency=afterok:$job2 7_visualization.sh | awk '{print $4}')
echo "Submitted 7_visualization.sh (JobID: $job3)"

job4=$(sbatch --dependency=afterok:$job3 8_advanced_analysis.sh | awk '{print $4}')
echo "Submitted 8_advanced_analysis.sh (JobID: $job4)"

job5=$(sbatch --dependency=afterok:$job4 9_gene_overlap_analysis.sh | awk '{print $4}')
echo "Submitted 9_gene_overlap_analysis.sh (JobID: $job5)"

job6=$(sbatch --dependency=afterok:$job5 10_go_enrichment.sh | awk '{print $4}')
echo "Submitted 10_go_enrichment.sh (JobID: $job6)"

echo "All jobs submitted. Pipeline will execute in sequence."
echo "Monitor progress with: squeue -u \$USER" 