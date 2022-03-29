# script to run with infinite time on west partition and go through all of the thermo calculations


import pandas as pd
import os
import sys
import time
import subprocess
import job_manager


# read in the species csv
print('Collecting Remaining Species')
try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'


species_csv = '/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/species_list.csv'
species_df = pd.read_csv(species_csv)


# for i in range(0, len(species_df)):
for i in [13]:
    species_index = species_df.i.values[i]
    species_smiles = species_df.SMILES.values[i]
    arkane_result = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'arkane', 'RMG_libraries', 'thermo.py')
    if os.path.exists(arkane_result):
        print(species_index, species_smiles, 'COMPLETE')
    else:
        print(species_index, species_smiles, 'RUNNING THERMO')
        start = time.time()
        end = time.time()

        # start a job that calls snakemake to run conformers
        conformer_cmd = f'snakemake -c1 species_thermo --config species_index={species_index}'
        print(f'Running {conformer_cmd}')
        cmd_pieces = conformer_cmd.split()
        self.proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
        print(self.proc)

        # hopefully, self.proc has the name of the next slurm job
        hotbit_job_num = 193829482


        # wait for the job to finish
        # Submitted batch job 24468402



        # wait for the job to 

        exit(0)
        # # start a job that calls snakemake to run rotors
        # conformer_cmd = f'snakemake -c1 run_rotors --config species_index={species_index}'

        # # start a job that calls snakemake to run arkane
        # conformer_cmd = f'snakemake -c1 run_arkane_thermo --config species_index={species_index}'

        # duration = end - start
        # print(f'COMPLETED IN {duration} SECONDS')