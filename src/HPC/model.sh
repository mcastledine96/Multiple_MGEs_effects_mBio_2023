#!/bin/bash -l
# Use bash and pickup a basic login environment.
#SBATCH --begin=now
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=defq                    # Specifies that the job will run on the default queue nodes.
#SBATCH --job-name=phage_fx       # A name for the job to be used when Job Monitoring
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL,
#SBATCH --mail-user=arn203@exeter.ac.uk   # E-Mail address of the user that needs to be notified.
#SBATCH --output=array_%A-%a.log                # Standard output and error log
#SBATCH --array=1-100                          # Array range
pwd; hostname; date

echo "Running a program on $SLURM_JOB_NODELIST"

module load Workspace/v1
export JULIA_NUM_THREADS=16
export mkdir ${WORKSPACE}/phage_plas/tmp
export TMPDIR=${WORKSPACE}/phage_plas/tmp
julia model_publication.jl $SLURM_ARRAY_TASK_ID --threads 16
julia model_si.jl $SLURM_ARRAY_TASK_ID --threads 16

seff.patched $SLURM_JOB_ID
