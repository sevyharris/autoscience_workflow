#!/bin/bash
#SBATCH --job-name=species
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --mincpus=32

# must pass the species index as the first argument
python /work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/species_thermo.py $1


