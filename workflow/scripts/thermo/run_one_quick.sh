#!/bin/bash
#SBATCH --job-name=one_thermo
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --partition=express,short,west
#SBATCH --time=00:20:00


cd "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/"
python "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/thermo/run_one.py" $1
