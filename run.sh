preset1cpu2gb="--partition=slims --cpus-per-task=1 --mem-per-cpu=2G"
sbatch $preset1cpu2gb -J experiment run_experiment_job.sh
