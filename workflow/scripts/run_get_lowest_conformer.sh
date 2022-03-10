#!/bin/bash
#SBATCH --job-name=lowest_conformer
#SBATCH --nodes=1

# must pass the species index as the first argument
python /work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/get_lowest_conformer.py $1

