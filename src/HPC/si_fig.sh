#!/bin/bash -l
# Use bash and pickup a basic login environment.
#SBATCH --begin=now
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=hmq                    # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=phage_fx       # A name for the job to be used when Job Monitoring
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL,
#SBATCH --mail-user=arn203@exeter.ac.uk   # E-Mail address of the user that needs to be notified.
#SBATCH --output=array_%A-%a.log                # Standard output and error log
#SBATCH --array=1-99                          # Array range
pwd; hostname; date

echo "Running a program on $SLURM_JOB_NODELIST"

module load Workspace/v1
export JULIA_NUM_THREADS=4
julia si_fig.jl $SLURM_ARRAY_TASK_ID --threads 4

seff.patched $SLURM_JOB_ID
