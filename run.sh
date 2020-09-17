make
mkdir -p results

preset1cpu2gb="--partition=slims --cpus-per-task=1 --mem-per-cpu=2200M"
sbatch $preset1cpu2gb -J experiment run_experiment_job.sh
