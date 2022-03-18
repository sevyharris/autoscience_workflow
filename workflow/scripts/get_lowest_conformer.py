# Script to generate species gaussian files/job

import os
import glob
import sys
import shutil

import numpy as np
import pandas as pd

import ase.io
import ase.io.gaussian
from ase.calculators.calculator import PropertyNotImplementedError

import rmgpy.species
import rmgpy.chemkin


import autotst.species
import autotst.calculator.gaussian


print_only = False  # option to just display the lowest energy conformer without copying
if len(sys.argv) > 2 and sys.argv[2] == 'print_only':
    print_only = True


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


# TODO check that the conformer calculation jobs are all done
# Get list of conformers
conformer_files = glob.glob(os.path.join(conformer_dir, 'conformer_*.log'))

energies = np.zeros(len(conformer_files))
for i, conformer_file in enumerate(conformer_files):
    try:
        with open(conformer_file, 'r') as f:
            atoms = ase.io.gaussian.read_gaussian_out(f)
            energy = atoms.get_potential_energy()
            energies[i] = energy
    except IndexError:
        pass
    except PropertyNotImplementedError:
        pass

lowest_idx = np.argmin(energies)
print(f'lowest energy one is {conformer_files[lowest_idx]}')
if print_only:
    exit(0)

# copy the lowest energy one into the rotor folder
os.makedirs(rotor_dir, exist_ok=True)
shutil.copy(conformer_files[lowest_idx], rotor_dir)

