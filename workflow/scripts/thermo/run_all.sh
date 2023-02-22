#!/bin/bash
#SBATCH --job-name=all_thermo
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --partition=short,west
#SBATCH --time=1-00:00:00
#SBATCH --mincpus=1
#SBATCH --exclude=c5003


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/thermo/run_all.py"
