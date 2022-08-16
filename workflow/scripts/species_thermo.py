# Script to generate species gaussian files/job

import os
import glob
import sys

import pandas as pd

import ase.io

import rmgpy.species
import rmgpy.chemkin

from hotbit import Hotbit

import autotst.species
from autotst.calculator.gaussian import Gaussian
import ase.calculators.lj

import job_manager


DFT_DIR = os.environ['DFT_DIR']
species_index = int(sys.argv[1])
print(f'Species index is {species_index}')


# Load the species from the official species list
scripts_dir = os.path.dirname(__file__)
# species_csv = os.path.join(scripts_dir, '..', '..', 'resources', 'species_list.csv')
species_csv = os.path.join(DFT_DIR, 'species_list.csv')
species_df = pd.read_csv(species_csv)

species_smiles = species_df.SMILES[species_index]
spec = autotst.species.Species([species_smiles])

print(f"loaded species {species_smiles}")
thermo_base_dir = os.path.join(DFT_DIR, 'thermo')
species_base_dir = os.path.join(thermo_base_dir, f'species_{species_index:04}')
os.makedirs(species_base_dir, exist_ok=True)


# generate conformers
try:
    spec.generate_conformers(ase_calculator=Hotbit())
    n_conformers = len(spec.conformers[species_smiles])
    print(f'{n_conformers} found with Hotbit')
except RuntimeError:
    # if hotbit fails, use built-in lennard jones
    print('Using built-in ase LennardJones calculator instead of Hotbit')
    spec.generate_conformers(ase_calculator=ase.calculators.lj.LennardJones())
    n_conformers = len(spec.conformers[species_smiles])
    print(f'{n_conformers} found with ase LennardJones calculator')

# do detailed calculation using Gaussian
conformer_dir = os.path.join(species_base_dir, 'conformers')
# write Gaussian input files
print("generating gaussian input files")
for i, cf in enumerate(spec.conformers[species_smiles]):
    gaussian = Gaussian(conformer=cf)
    calc = gaussian.get_conformer_calc()
    calc.label = f'conformer_{i:04}'
    calc.directory = conformer_dir
    calc.parameters.pop('scratch')
    calc.parameters.pop('multiplicity')
    calc.parameters['mult'] = cf.rmg_molecule.multiplicity
    calc.chk = f'conformer_{i:04}.chk'
    calc.write_input(cf.ase_molecule)


# Make slurm script
# Make a file to run Gaussian
slurm_run_file = os.path.join(conformer_dir, 'run.sh')
slurm_settings = {
    '--job-name': f'g16_cf_{species_index}',
    '--error': 'error.log',
    '--nodes': 1,
    '--partition': 'west,short',
    '--exclude': 'c5003',
    '--mem': '20Gb',
    '--time': '24:00:00',
    '--cpus-per-task': 16,
    '--array': f'0-{n_conformers - 1}%30',
}

slurm_file_writer = job_manager.SlurmJobFile(
    full_path=slurm_run_file,
)
slurm_file_writer.settings = slurm_settings

slurm_file_writer.content = [
    'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
    'mkdir -p $GAUSS_SCRDIR\n',
    'module load gaussian/g16\n',
    'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

    'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
    'fname="conformer_${RUN_i}.com"\n\n',

    'g16 $fname\n',
]
slurm_file_writer.write_file()

# submit the job
start_dir = os.getcwd()
os.chdir(conformer_dir)
gaussian_conformers_job = job_manager.SlurmJob()
slurm_cmd = f"sbatch {slurm_run_file}"
gaussian_conformers_job.submit(slurm_cmd)
os.chdir(start_dir)
