#!/bin/bash
#SBATCH --job-name=kinetics_job
#SBATCH --error=kinetics_error.log
#SBATCH --mincpus=1
#SBATCH --exclude=c5003
#SBATCH --partition=west,short,express
#SBATCH --time=00:30:00


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/kinetics/run_one.py" $1
