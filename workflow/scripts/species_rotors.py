import os
import sys
import glob

import pandas as pd
import ase.io.gaussian

import autotst.species
import autotst.calculator.gaussian

import job_manager

import rotor_scan

# Read in the species
DFT_DIR = os.environ['DFT_DIR']
species_index = int(sys.argv[1])
print(f'Species index is {species_index}')


# Load the species from the official species list
scripts_dir = os.path.dirname(__file__)
species_csv = os.path.join(scripts_dir, '..', '..', 'resources', 'species_list.csv')
species_df = pd.read_csv(species_csv)
species_smiles = species_df.SMILES[species_index]

print(f"loaded species {species_smiles}")
thermo_base_dir = os.path.join(DFT_DIR, 'thermo')
species_base_dir = os.path.join(thermo_base_dir, f'species_{species_index:04}')
conformer_dir = os.path.join(species_base_dir, 'conformers')
rotor_dir = os.path.join(species_base_dir, 'rotors')

# confirm there's only one conformer file
conformer_files = glob.glob(os.path.join(rotor_dir, 'conformer_*.log'))
if len(conformer_files) > 1:
    print(f"Warning: more than one lowest energy conformer. Using {conformer_files[0]}")

with open(conformer_files[0], 'r') as f:
    atoms = ase.io.gaussian.read_gaussian_out(f)

# make a conformer object again
new_cf = autotst.species.Conformer(smiles=species_smiles)
new_cf._ase_molecule = atoms
new_cf.update_coords_from(mol_type="ase")

# get the rotors
torsions = new_cf.get_torsions()
n_rotors = len(torsions)

print("generating gaussian input files")
# gaussian = autotst.calculator.gaussian.Gaussian(conformer=new_cf)
for i, torsion in enumerate(new_cf.torsions):
    # print(torsion)
    # calc = gaussian.get_rotor_calc(torsion_index=i)
    # calc.label = f'rotor_{i:04}'
    # calc.directory = rotor_dir
    # calc.parameters.pop('scratch')
    # calc.parameters.pop('multiplicity')
    # calc.parameters['mult'] = new_cf.rmg_molecule.multiplicity
    # calc.write_input(new_cf.ase_molecule)

    fname = os.path.join(rotor_dir, f'rotor_{i:04}.com')
    rotor_scan.write_scan_file(fname, new_cf, i)


# Make a slurm script to run all rotors
slurm_run_file = os.path.join(rotor_dir, 'run_rotor_calcs.sh')
slurm_settings = {
    '--job-name': f'g16_rotors_{species_index}',
    '--error': 'error.log',
    '--nodes': 1,
    '--partition': 'west,short',
    '--mem': '20Gb',
    '--time': '24:00:00',
    '--cpus-per-task': 16,
    '--array': f'0-{n_rotors - 1}%40',
}

slurm_file_writer = job_manager.SlurmJobFile(
    full_path=slurm_run_file,
)
slurm_file_writer.settings = slurm_settings
slurm_file_writer.content = [
    # TODO figure out if I really need to include the GAUSS_SCRDIR
    'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
    'mkdir -p $GAUSS_SCRDIR\n',
    'module load gaussian/g16\n',
    'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

    'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
    'fname="rotor_${RUN_i}.com"\n\n',

    'g16 $fname\n',
]
slurm_file_writer.write_file()

# submit the job
start_dir = os.getcwd()
os.chdir(rotor_dir)
gaussian_rotors_job = job_manager.SlurmJob()
slurm_cmd = f"sbatch {slurm_run_file}"
gaussian_rotors_job.submit(slurm_cmd)
os.chdir(start_dir)


