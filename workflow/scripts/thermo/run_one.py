# Script for running a thermo job for one species
import sys
import job


species_index = int(sys.argv[1])
species_smiles = job.index2smiles(species_index)
print(species_smiles)

if job.arkane_complete(species_index):
    print("COMPLETE")
    exit(0)

print("RUNNING CALCULATION")
job.run_conformers_job(species_index)
job.run_rotors_job(species_index)
job.run_arkane_job(species_index)
