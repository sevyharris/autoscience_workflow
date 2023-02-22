#!/bin/bash
#SBATCH --job-name=kinetics_job
#SBATCH --partition=west,short
#SBATCH --time=24:00:00


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/kinetics/run_one.py" $1
