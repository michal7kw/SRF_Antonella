#!/bin/bash
#SBATCH --parsable
#SBATCH --partition={resources.slurm_partition}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={resources.mem_mb}M
#SBATCH --time={resources.runtime}
#SBATCH --output=logs/cluster_logs/%j.out
#SBATCH --error=logs/cluster_logs/%j.err

# Execute the command
{exec_job} 