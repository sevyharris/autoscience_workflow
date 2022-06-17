# script to copy lowest conformer to the rotors folder
import os
import sys
import job
import shutil


species_index = int(sys.argv[1])
best_conformer_file = job.get_lowest_conformer(species_index)
if best_conformer_file is None:
    raise ValueError('No valid conformers!')

rotor_dir = os.path.join(job.DFT_DIR, 'thermo', f'species_{species_index:04}', 'rotors')
os.makedirs(rotor_dir, exist_ok=True)
shutil.copy(best_conformer_file, rotor_dir)
