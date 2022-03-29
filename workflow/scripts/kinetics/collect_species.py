# Script to figure out which species are needed to run kinetics
import os
import sys

import pandas as pd

# import ase.io

# import rmgpy.species
# import rmgpy.chemkin


try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    DFT_DIR = "/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft"
reaction_index = int(sys.argv[1])
print(f'Preparing reaction {reaction_index}')

# Load the species from the official species list

scripts_dir = os.path.dirname(os.path.dirname(__file__))
print(scripts_dir)
reaction_csv = os.path.join(scripts_dir, '..', '..', 'resources', 'reaction_list.csv')
reaction_df = pd.read_csv(reaction_csv)
species_csv = os.path.join(scripts_dir, '..', '..', 'resources', 'species_list.csv')
species_df = pd.read_csv(species_csv)

reaction_smiles = reaction_df.SMILES[reaction_index]
print(f'Collecting relevant species for reaction {reaction_index}: {reaction_smiles}')

reactants = reaction_smiles.split('_')[0].split('+')
products = reaction_smiles.split('_')[1].split('+')
species_list = reactants + products


incomplete_species = []
print("Required Species:")
for spec in species_list:
    spec_index = species_df.index[species_df['SMILES'] == spec]
    row = species_df.loc[spec_index]
    
    print(row.i.values[0], row.SMILES.values[0])
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{row.i.values[0]:04}')
    if not os.path.exists(os.path.join(species_dir, 'arkane')):
        incomplete_species.append([row.i.values[0], row.SMILES.values[0]])

print('Missing Species:')
for spec in incomplete_species:
    print(spec)
