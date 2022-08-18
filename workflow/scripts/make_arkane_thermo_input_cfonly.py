# No rotors, don't need to use autotst

import os
import sys
import glob
import shutil

import cclib.io
import pandas as pd
import ase.io.gaussian

import job_manager

import rmgpy.species


# Read in the species
DFT_DIR = os.environ['DFT_DIR']
species_index = int(sys.argv[1])
print(f'Species index is {species_index}')


# Load the species from the official species list
scripts_dir = os.path.dirname(__file__)
species_csv = os.path.join(DFT_DIR, 'species_list.csv')
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
elif len(conformer_files) == 0:
    print(f'Conformer files not found. Did you remember to run the rotors?')

# copy to the arkane folder
arkane_dir = os.path.join(species_base_dir, 'arkane')
os.makedirs(arkane_dir, exist_ok=True)
shutil.copy(conformer_files[0], arkane_dir)
conformer_file = os.path.join(arkane_dir, os.path.basename(conformer_files[0]))


with open(conformer_file, 'r') as f:
    atoms = ase.io.gaussian.read_gaussian_out(f)

# get around making a conformer object from the SMILES
rmg_molecule = rmgpy.species.Species(smiles=species_smiles)


def write_conformer_file(rmg_molecule, gauss_log, arkane_dir, include_rotors=False):
    # assume rotor and conformer logs have already been copied into the arkane directory
    label = rmg_molecule.smiles
    species_name = os.path.basename(gauss_log[:-4])
    # parser = cclib.io.ccread(gauss_log)
    # symbol_dict = {
    #     35: "Br",
    #     17: "Cl",
    #     9:  "F",
    #     8:  "O",
    #     7:  "N",
    #     6:  "C",
    #     1:  "H",
    #     18: "Ar",
    #     2: "He",
    #     10: "Ne",
    # }

    # atoms = []

    # for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
    #     atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))

    # conformer._ase_molecule = ase.Atoms(atoms)
    # conformer.update_coords_from("ase")
    # mol = conformer.rmg_molecule
    output = ['#!/usr/bin/env python',
                  '# -*- coding: utf-8 -*-', ]

    output += ["",
               f"spinMultiplicity = {rmg_molecule.multiplicity}",
               ""]
    model_chemistry = 'M06-2X/cc-pVTZ'

    # use relative path for easy transfer -- assume we will copy the log files into the Arkane folder
    gauss_log_relative = os.path.basename(gauss_log)
    output += ["energy = {", f"    '{model_chemistry}': Log('{gauss_log_relative}'),", "}", ""]  # fix this

    output += [f"geometry = Log('{gauss_log_relative}')", ""]
    output += [
        f"frequencies = Log('{gauss_log_relative}')", ""]

    input_string = ""

    for t in output:
        input_string += t + "\n"

    with open(os.path.join(arkane_dir, species_name + '.py'), "w") as f:
        f.write(input_string)
    return True


# write the Arkane conformer file
write_conformer_file(rmg_molecule.molecule[0], conformer_file, arkane_dir, include_rotors=False)

# write the Arkane input file
input_file = os.path.join(arkane_dir, 'input.py')
formula = rmg_molecule.molecule[0].get_formula()
lines = [
    '#!/usr/bin/env python\n\n',
    f'modelChemistry = "M06-2X/cc-pVTZ"\n',
    'useHinderedRotors = True\n',
    'useBondCorrections = False\n\n',

    'frequencyScaleFactor = 0.982\n',

    f"species('{formula}', '{os.path.basename(conformer_file[:-4])}.py', structure=SMILES('{rmg_molecule.smiles}'))\n\n",

    f"thermo('{formula}', 'NASA')\n",
]
with open(input_file, 'w') as f:
    f.writelines(lines)


# copy a run script into the arkane directory
run_script = os.path.join(arkane_dir, 'run_arkane.sh')
with open(run_script, 'w') as f:
    # Run on express
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --partition=express,short,west\n')
    f.write('#SBATCH --time=00:20:00\n\n')
    f.write('python ~/rmg/RMG-Py/Arkane.py input.py\n\n')
