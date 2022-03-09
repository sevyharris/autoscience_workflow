# Script to generate species gaussian files/job

import os
import glob
import sys

import ase.io

import rmgpy.species
import rmgpy.chemkin

from hotbit import Hotbit

import autotst.species
from autotst.calculator.gaussian import Gaussian

import job_manager



# Load the model:
# load the model
chemkin_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/chem_annotated.inp"
dictionary_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/species_dictionary.txt"
transport_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/tran.dat"
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin_path,
    dictionary_path=dictionary_path,
    transport_path=transport_path
)
print (f"Loaded model with {len(species_list)} species and {len(reaction_list)} reactions")



# TODO connect ranked list to this code. for now, just change up the species index
# TODO connect the ranked list to this code
species_index = 80  # 3 corresponds to n-heptane, the top species
species_index = 79  # 3 corresponds to n-heptane, the top species
species_index = 71  # incomplete
species_index = 68  # incomplete

species_index = int(sys.argv[1])

spec_rmg = species_list[species_index]
spec = autotst.species.Species([spec_rmg.smiles])

print(f"loaded species {spec_rmg}")
thermo_base_dir = '/work/westgroup/harris.se/autoscience/dft/thermo'
species_base_dir = os.path.join(thermo_base_dir, f'species_{species_index:04}')
os.makedirs(species_base_dir, exist_ok=True)


# generate conformers
spec.generate_conformers(ase_calculator=Hotbit())
n_conformers = len(spec.conformers[spec_rmg.smiles])
print(f'{n_conformers} found with Hotbit')


# do detailed calculation using Gaussian
conformer_dir = os.path.join(species_base_dir, 'conformers')
# write Gaussian input files
print("generating gaussian input files")
for i, cf in enumerate(spec.conformers[spec_rmg.smiles]):
    gaussian = Gaussian(conformer=cf)
    calc = gaussian.get_conformer_calc()
    calc.label = f'conformer_{i:04}'
    calc.directory = conformer_dir
    calc.parameters.pop('scratch')
    calc.parameters.pop('multiplicity')
    calc.parameters['mult'] = cf.rmg_molecule.multiplicity
    calc.write_input(cf.ase_molecule)


# Make slurm script
# Make a file to run Gaussian
slurm_run_file = os.path.join(conformer_dir, 'run.sh')
slurm_settings = {
    '--job-name': f'g16_cf_{species_index}',
    '--error': 'error.log',
    '--output': 'output.log',
    '--nodes': 1,
    '--partition': 'west,short',
    '--exclude': 'c5003',
    '--mem': '20Gb',
    '--time': '24:00:00',
    '--cpus-per-task': 16,
    '--array': f'0-{n_conformers - 1}%40',
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





