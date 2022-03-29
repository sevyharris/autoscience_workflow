#!/bin/bash
#SBATCH --job-name=ghost_g16_thermo
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --partition=west
#SBATCH --time=7-00:00:00
#SBATCH --mincpus=1
#SBATCH --exclude=c5003


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/run_all_thermo.py"
