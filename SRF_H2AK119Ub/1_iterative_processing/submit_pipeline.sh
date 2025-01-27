#!/bin/bash

# Submit first job and capture its ID
job1_id=$(sbatch --parsable 5_differential_binding.sh)
echo "Submitted differential binding analysis job: $job1_id"

# Submit annotation job dependent on first job
job2_id=$(sbatch --parsable --dependency=afterok:${job1_id} 6_annotation_and_enrichment.sh)
echo "Submitted annotation job: $job2_id"

# Submit visualization job dependent on annotation job
job3_id=$(sbatch --parsable --dependency=afterok:${job2_id} 7_visualization.sh)
echo "Submitted visualization job: $job3_id"

# Submit advanced analysis job dependent on visualization job
# job4_id=$(sbatch --parsable --dependency=afterok:${job3_id} 9_advanced_analysis.sh)
# echo "Submitted advanced analysis job: $job4_id"

echo "All jobs submitted successfully"