#!/bin/bash

# Submit first job and capture its ID
job1_id=$(sbatch --parsable ./sh/2a_quality_control.sh)
echo "Submitted quality control job: $job1_id"

# Submit annotation job dependent on first job
job2_id=$(sbatch --parsable --dependency=afterok:${job1_id} ./sh/2b_multiqc.sh)
echo "Submitted multiqc job: $job2_id"

# Submit visualization job dependent on annotation job
job3_id=$(sbatch --parsable --dependency=afterok:${job2_id} ./sh/3_alignment.sh)
echo "Submitted alignment job: $job3_id"

# Submit advanced analysis job dependent on visualization job
job4_id=$(sbatch --parsable --dependency=afterok:${job3_id} ./sh/4_peak_calling2_improved.sh)
echo "Submitted peak calling job: $job4_id"

# Submit advanced analysis job dependent on visualization job
job5_id=$(sbatch --parsable --dependency=afterok:${job4_id} ./sh/4b_bigwig.sh)
echo "Submitted bigwig job: $job5_id"

echo "All jobs submitted successfully"