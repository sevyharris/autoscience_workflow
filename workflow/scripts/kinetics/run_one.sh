#!/bin/bash
#SBATCH --job-name=kinetics_job
#SBATCH --error=kinetics_error.log
#SBATCH --mincpus=1
#SBATCH --exclude=c5003
#SBATCH --partition=west,short
#SBATCH --time=1-00:00:00


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/kinetics/run_one.py" $1
