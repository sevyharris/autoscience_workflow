# Script for running a thermo job for one species
import sys
import job


logfile = 'all_logs.txt'
n_species = job.get_num_species()
print(n_species)
# skip_indices = [1, 2, 3, 4, 9, 10, 42, 45, 46, 47]
skip_indices = [9, 10]
for species_index in range(110, 180):
    species_smiles = job.index2smiles(species_index)
    if species_index in skip_indices:
        with open(logfile, 'a') as f:
            f.write(f'Skipping species {species_index}: {species_smiles}' + '\n')
        print(f'Skipping species {species_index}: {species_smiles}')
        continue

    if job.arkane_complete(species_index):
        with open(logfile, 'a') as f:
            f.write(f"SPECIES {species_index}: {species_smiles} already COMPLETE" + '\n')
        print(f"SPECIES {species_index}: {species_smiles} already COMPLETE")
        continue

    with open(logfile, 'a') as f:
        f.write(f"Running Calculation for Species {species_index}: {species_smiles}" + '\n')
    print(f"Running Calculation for Species {species_index}: {species_smiles}")
    job.run_conformers_job(species_index)
    job.run_rotors_job(species_index)
    job.run_arkane_job(species_index)
    with open(logfile, 'a') as f:
        f.write(f"SPECIES {species_index}: {species_smiles} COMPLETE" + '\n')
    print(f"SPECIES {species_index}: {species_smiles} COMPLETE")
