#!/bin/bash
#SBATCH --job-name=one_thermo
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --partition=west
#SBATCH --time=7-00:00:00
#SBATCH --mincpus=1
#SBATCH --exclude=c5003


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/thermo/run_one.py" $1
