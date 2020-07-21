preset1cpu2gb="--partition=general --cpus-per-task=1 --mem-per-cpu=4000M"
sbatch $preset1cpu2gb -J experiment run_experiment_job.sh
