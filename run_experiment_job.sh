#!/bin/bash -e
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=franciscojacb@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=1-100   ## Separate into different tasks

fname=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
./bin/experiment $SLURM_ARRAY_TASK_ID > results/"$fname"

