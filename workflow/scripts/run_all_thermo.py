# script to run with infinite time on west partition and go through all of the thermo calculations


import pandas as pd
import os
import sys
import glob
import time
import subprocess
import job_manager


logfile = 'ghost.log'
# read in the species csv
with open(logfile, 'a') as f:
    f.write('Collecting Remaining Species\n')
print('Collecting Remaining Species')

try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'


workflow_dir = '/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/'
species_csv = '/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/species_list.csv'
species_df = pd.read_csv(species_csv)


for i in range(0, len(species_df)):
# for i in [13]:
    species_index = species_df.i.values[i]
    species_smiles = species_df.SMILES.values[i]
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    arkane_result = os.path.join(species_dir, 'arkane', 'RMG_libraries', 'thermo.py')
    if os.path.exists(arkane_result):
        with open(logfile, 'a') as f:
            f.write(species_index, species_smiles, 'COMPLETE\n')
        print(species_index, species_smiles, 'COMPLETE')
    else:
        print(species_index, species_smiles, 'RUNNING THERMO')
        start = time.time()
        end = time.time()

        # start a job that calls snakemake to run conformers
        os.chdir(workflow_dir)
        conformer_cmd = f'snakemake -c1 species_thermo --config species_index={species_index}'
        with open(logfile, 'a') as f:
            f.write(f'Running {conformer_cmd}' + '\n')
        print(f'Running {conformer_cmd}')
        cmd_pieces = conformer_cmd.split()
        proc = subprocess.Popen(cmd_pieces)
        print(proc)


        # wait 10 minutes for the hotbit job to start/finish
        # try to read the slurm file in 
        time.sleep(600)
        g16_job_number = ''
        hotbit_slurm = glob.glob(os.path.join(species_dir, 'slurm-*'))
        if len(hotbit_slurm) == 0:
            print('Hotbit slurm file not found')
            exit(3)
        hotbit_complete = False
        while not hotbit_complete:
            with open(hotbit_slurm[0], 'r') as f:
                lines = f.readlines()
                for line in line:
                    if 'Submitted batch job' in line:
                        hotbit_complete = True
                        g16_job_number = line.split()[-1]
                        break
            time.sleep(600)
        print('Hotbit conformer screening complete')
        with open(logfile, 'a') as f:
            f.write('Hotbit conformer screening complete\n')


        # wait 10 minutes for the conformer jobs to finish
        gaussian_job = job_manager.SlurmJob()
        gaussian_job.job_id = g16_job_number
        print(f'Waiting on job {gaussian_job}')
        with open(logfile, 'a') as f:
            f.write(f'Waiting on job {g16_job_number}' + '\n')
        gaussian_job.wait_all(check_interval=3600)
        print('Gaussian jobs complete')
        with open(logfile, 'a') as f:
            f.write('Gaussian jobs complete\n')


        # start a job that calls snakemake to run rotors
        rotor_cmd = f'snakemake -c1 run_rotors --config species_index={species_index}'
        print(f'Running {rotor_cmd}')
        cmd_pieces = rotor_cmd.split()
        proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
        print(proc)

        # wait 10 minutes for the rotor gaussian job to begin
        time.sleep(600)
        g16_job_number = ''
        
        skip_rotors = False
        if os.path.exists(os.path.join(species_dir, 'rotors', 'NO_ROTORS.txt')):
            skip_rotors = True

        if not skip_rotors:
            rotor_slurm_files = glob.glob(os.path.join(species_dir, 'rotors', 'slurm-*'))
            if len(rotor_slurm_files) == 0:
                print('Rotor slurm file not found')
                exit(3)
            
            rotor_slurm_file = os.path.basename(rotor_slurm_files[0])
            rotor_slurm_id = rotor_slurm_file[6:14]
            rotor_job = job_manager.SlurmJob()
            rotor_job.job_id = rotor_slurm_id
            print(f'Waiting on job {rotor_slurm_id}')
            with open(logfile, 'a') as f:
                f.write(f'Waiting on job {rotor_slurm_id}' + '\n')
            rotor_job.wait_all(check_interval=3600)
            print('Rotor jobs complete')
            with open(logfile, 'a') as f:
                f.write('Rotor jobs complete\n')


        # start a job that calls snakemake to run arkane
        arkane_cmd = f'snakemake -c1 run_arkane_thermo --config species_index={species_index}'
        print(f'Running {arkane_cmd}')
        cmd_pieces = arkane_cmd.split()
        proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
        print(proc)

        
        
        # wait 10 minutes for Arkane start/finish
        # try to read the slurm file in
        print('Waiting for arkane job')
        with open(logfile, 'a') as f:
            f.write('Waiting for arkane job\n')
        while not os.path.exists(arkane_result):
            time.sleep(600)
        print('Arkane complete')
        with open(logfile, 'a') as f:
            f.write('Arkane complete\n')

        duration = end - start
        print(f'COMPLETED IN {duration} SECONDS')
        with open(logfile, 'a') as f:
            f.write(f'COMPLETED IN {duration} SECONDS' + '\n')
