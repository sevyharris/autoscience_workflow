import os
import glob
import shutil

import arkane.ess.gaussian
import arkane.exceptions
import autotst.reaction
import rmgpy.chemkin

# from hotbit import Hotbit

import sys
import kineticfun

sys.path.append('/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/thermo/')
sys.path.append('/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/kinetics/')
import job


def get_reaction_label(rmg_reaction):
    label = ''
    for reactant in rmg_reaction.reactants:
        label += get_sp_name(reactant.smiles) + ' + '
    label = label[:-2]
    label += '<=> '
    for product in rmg_reaction.products:
        label += get_sp_name(product.smiles) + ' + '
    label = label[:-3]
    return label


# Read in the reaction

reaction_index = int(sys.argv[1])
reaction_smiles = kineticfun.reaction_index2smiles(reaction_index)
print(f'Reaction index is {reaction_index}')

# DFT_DIR = os.environ['DFT_DIR']
DFT_DIR = "/work/westgroup/harris.se/autoscience/autoscience/butane/dft"
reaction_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
overall_dir = os.path.join(reaction_dir, 'overall')
arkane_dir = os.path.join(reaction_dir, 'arkane')
os.makedirs(arkane_dir, exist_ok=True)

species_dict_file = "/work/westgroup/harris.se/autoscience/autoscience/butane/models/rmg_model/species_dictionary.txt"
species_dict = rmgpy.chemkin.load_species_dictionary(species_dict_file)


def get_sp_name(smiles):
    if smiles == '[CH2]C=CC':  # manually change to resonance structures included in model
        smiles = 'C=C[CH]C'
    elif smiles == '[CH2][CH]C=C':
        smiles = '[CH2]C=C[CH2]'
    for entry in species_dict.keys():
        if species_dict[entry].smiles == smiles:
            return str(species_dict[entry])
    # need to look for isomorphism
    print(f'Failed to get species name for {smiles}')


direction = 'forward'
reaction = autotst.reaction.Reaction(label=reaction_smiles)
reaction.ts[direction][0].get_molecules()

# check whether autotst renamed a species to one of its resonance structures to get a match

# reaction.generate_conformers(ase_calculator=Hotbit())  # probably have to add this back in for multiple TS handling


# pick the lowest energy transition state:
TS_logs = glob.glob(os.path.join(overall_dir, f'fwd_ts_*.log'))
TS_log = ''
lowest_energy = 0
for logfile in TS_logs:
    try:
        g_reader = arkane.ess.gaussian.GaussianLog(logfile)
        energy = g_reader.load_energy()
        if energy < lowest_energy:
            lowest_energy = energy
            TS_log = logfile
    except arkane.exceptions.LogError:
        print(f'skipping bad logfile {logfile}')
        continue

# ----------------------------------------------------------------- #
# write the input file
model_chemistry = 'M06-2X/cc-pVTZ'

lines = [
    f'modelChemistry = "{model_chemistry}"\n',
    'useHinderedRotors = False\n',
    'useBondCorrections = False\n\n',
]

completed_species = []
for reactant in reaction.rmg_reaction.reactants + reaction.rmg_reaction.products:
    # check for duplicates
    duplicate = False
    for sp in completed_species:
        if reactant.is_isomorphic(sp):
            duplicate = True
    if duplicate:
        continue

    species_smiles = reactant.smiles
    if species_smiles == '[CH2]C=CC':  # TODO clean up this fix where we manually switch back to other resonance structure
        species_smiles = 'C=C[CH]C'
    elif species_smiles == '[CH2][CH]C=C':
        species_smiles = '[CH2]C=C[CH2]'
    species_name = get_sp_name(species_smiles)
    species_index = job.smiles2index(species_smiles)
    species_arkane_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'arkane')

    species_file = os.path.join(f'species_{species_index:04}', os.path.basename(glob.glob(os.path.join(species_arkane_dir, 'conformer_*.py'))[0]))

    try:
        shutil.copytree(species_arkane_dir, os.path.join(arkane_dir, f'species_{species_index:04}'))
    except FileExistsError:
        pass

    # TODO - copy the species into the destination so the arkane calculation can be copied and redone elsewhere
    lines.append(f'species("{species_name}", "{species_file}", structure=SMILES("{species_smiles}"))\n')
    lines.append(f'thermo("{species_name}", "NASA")\n\n')

    completed_species.append(reactant)


lines.append('\n')

TS_name = 'TS'
TS_file = 'TS.py'
TS_arkane_path = os.path.join(arkane_dir, TS_file)
shutil.copy(TS_log, arkane_dir)

lines.append(f'transitionState("{TS_name}", "{TS_file}")\n')

reaction_label = get_reaction_label(reaction.rmg_reaction)
reactants = [get_sp_name(reactant.smiles) for reactant in reaction.rmg_reaction.reactants]
products = [get_sp_name(product.smiles) for product in reaction.rmg_reaction.products]
lines.append(f'reaction(\n')
lines.append(f'    label = "{reaction_label}",\n')
lines.append(f'    reactants = {reactants},\n')
lines.append(f'    products = {products},\n')
lines.append(f'    transitionState = "{TS_name}",\n')
lines.append(f'#    tunneling = "Eckart",\n')
lines.append(f')\n\n')

lines.append(f'statmech("{TS_name}")\n')
lines.append(f'kinetics("{reaction_label}")\n\n')

# write the TS file
ts_lines = [
    'energy = {"' + f'{model_chemistry}": Log("{os.path.basename(TS_log)}")' + '}\n\n',
    'geometry = Log("' + f'{os.path.basename(TS_log)}")' + '\n\n',
    'frequencies = Log("' + f'{os.path.basename(TS_log)}")' + '\n\n',
]
with open(TS_arkane_path, 'w') as g:
    g.writelines(ts_lines)


arkane_input_file = os.path.join(arkane_dir, 'input.py')
with open(arkane_input_file, 'w') as f:
    f.writelines(lines)
# ----------------------------------------------------------------- #

# make the slurm script to run arkane
run_script = os.path.join(arkane_dir, 'run_arkane.sh')
with open(run_script, 'w') as f:
    # Run on express
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --partition=express,short,west\n')
    f.write('#SBATCH --time=00:20:00\n\n')
    f.write('python ~/rmg/RMG-Py/Arkane.py input.py\n\n')
