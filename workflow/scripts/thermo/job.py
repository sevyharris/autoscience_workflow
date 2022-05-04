# Functions for running a thermo job using this workflow


import pandas as pd
import os
import sys
import glob
import time
import subprocess
import job_manager


try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'


def index2smiles(species_index):
    """Function to return species smiles given a species index
    looks up the results in the species_list.csv
    """
    species_csv = os.path.join(DFT_DIR, '..', '..', 'resources', 'species_list.csv')
    species_df = pd.read_csv(species_csv)
    species_smiles = species_df.SMILES.values[species_index]
    return species_smiles


def arkane_complete(species_index):
    """Function to check whether the arkane job is complete for a species
    Expects to find the following directory structure:
    DFT_DIR/thermo/species_XXXX/arkane/RMG_libraries/thermo.py
    Returns True if complete, False otherwise
    """
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    arkane_result = os.path.join(species_dir, 'arkane', 'RMG_libraries', 'thermo.py')
    return os.path.exists(arkane_result)


def incomplete_conformers(species_index):
    """Returns a list of indices of incomplete conformers that need to be rerun"""
    conformer_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'conformers')

    # Get #conformers from the array job script
    slurm_array_file = os.path.join(conformer_dir, 'run.sh')
    with open(slurm_array_file, 'r') as f:
        for line in f:
            if 'SBATCH --array=' in line:
                token = line.split('-')[-1]
                n_conformers = 1 + int(token.split('%')[0])
                break

    incomplete_cfs = []
    for cf_index in range(0, n_conformers):
        conformer_file = os.path.join(conformer_dir, f'conformer_{cf_index:04}.log')
        if not os.path.exists(conformer_file):
            incomplete_cfs.append(cf_index)
            continue
        with open(conformer_file, 'rb') as f:
            try:
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            last_line = f.readline().decode()
            if 'Normal termination' not in last_line:
                incomplete_cfs.append(cf_index)
    return incomplete_cfs


def conformers_complete(species_index):
    """Function to check whether all of the Gaussian conformer jobs have finished running.
    Looks at the run.sh script to find the highest conformer index, then searches each .log file
    for Normal termination
    """
    if incomplete_conformers(species_index):
        return False
    return True


def run_conformers_job(species_index, logfile=None):
    """Function to call snakemake rule to run conformers
    This function waits until all SLURM jobs are done, so it could take days
    """
    conformers_complete_file = os.path.join(species_dir, 'conformers', 'CONFORMERS_COMPLETE')  # file to signal completion
    if os.path.exists(conformers_complete_file):
        print('Conformers already ran')
        return True

    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    workflow_dir = os.path.join(DFT_DIR, '..', '..', 'workflow')

    # start a job that calls snakemake to run conformers
    os.chdir(workflow_dir)
    conformer_cmd = f'snakemake -c1 species_thermo --config species_index={species_index}'
    print(f'Running {conformer_cmd}')
    cmd_pieces = conformer_cmd.split()
    proc = subprocess.Popen(cmd_pieces)
    print(proc)

    # RUN HOTBIT
    time.sleep(300)
    g16_job_number = ''
    # look for the hotbit slurm file
    hotbit_slurm = glob.glob(os.path.join(species_dir, 'slurm-*'))
    if len(hotbit_slurm) == 0:
        print('Hotbit slurm file not found. Hotbit did not start.')
        exit(3)
    hotbit_complete = False
    while not hotbit_complete:
        with open(hotbit_slurm[0], 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'Submitted batch job' in line:
                    hotbit_complete = True
                    g16_job_number = line.split()[-1]
                    break
        time.sleep(300)  # This wait is to make sure the job is on the SLURM queue
    print('Hotbit conformer screening complete')
    with open(logfile, 'a') as f:
        f.write('Hotbit conformer screening complete\n')

    # wait 10 minutes for the conformer jobs to finish
    gaussian_job = job_manager.SlurmJob()
    gaussian_job.job_id = g16_job_number
    print(f'Waiting on job {gaussian_job}')
    with open(logfile, 'a') as f:
        f.write(f'Waiting on job {g16_job_number}' + '\n')
    gaussian_job.wait_all(check_interval=600)
    print('Gaussian jobs complete')
    with open(logfile, 'a') as f:
        f.write('Gaussian jobs complete\n')
    with open(conformers_complete_file, 'w') as f:
        f.write('CONFORMERS COMPLETE')


#     # start a job that calls snakemake to run rotors
#     rotor_cmd = f'snakemake -c1 run_rotors --config species_index={species_index}'
#     print(f'Running {rotor_cmd}')
#     cmd_pieces = rotor_cmd.split()
#     proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
#     print(proc)

#     # wait 5 minutes for the rotor gaussian job to begin
#     time.sleep(300)
#     g16_job_number = ''
    
#     skip_rotors = False
#     if os.path.exists(os.path.join(species_dir, 'rotors', 'NO_ROTORS.txt')):
#         skip_rotors = True

#     if not skip_rotors:
#         rotor_slurm_files = glob.glob(os.path.join(species_dir, 'rotors', 'slurm-*'))
#         if len(rotor_slurm_files) == 0:
#             print('Rotor slurm file not found')
#             exit(3)
        
#         rotor_slurm_file = os.path.basename(rotor_slurm_files[0])
#         rotor_slurm_id = rotor_slurm_file[6:14]
#         rotor_job = job_manager.SlurmJob()
#         rotor_job.job_id = rotor_slurm_id
#         print(f'Waiting on job {rotor_slurm_id}')
#         with open(logfile, 'a') as f:
#             f.write(f'Waiting on job {rotor_slurm_id}' + '\n')
#         rotor_job.wait_all(check_interval=600)
#         print('Rotor jobs complete')
#         with open(logfile, 'a') as f:
#             f.write('Rotor jobs complete\n')


#     # start a job that calls snakemake to run arkane
#     arkane_cmd = f'snakemake -c1 run_arkane_thermo --config species_index={species_index}'
#     print(f'Running {arkane_cmd}')
#     cmd_pieces = arkane_cmd.split()
#     proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
#     print(proc)

    
    
#     # wait 10 minutes for Arkane start/finish
#     # try to read the slurm file in
#     print('Waiting for arkane job')
#     with open(logfile, 'a') as f:
#         f.write('Waiting for arkane job\n')
#     while not os.path.exists(arkane_result):
#         time.sleep(300)
#     print('Arkane complete')
#     with open(logfile, 'a') as f:
#         f.write('Arkane complete\n')

#     end = time.time()
#     duration = end - start
#     print(f'COMPLETED {species_smiles} IN {duration} SECONDS')
#     with open(logfile, 'a') as f:
#         f.write(f'COMPLETED {species_smiles} IN {duration} SECONDS' + '\n')
