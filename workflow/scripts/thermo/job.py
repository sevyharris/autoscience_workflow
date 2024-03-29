# Functions for running a thermo job using this workflow
import pandas as pd
import os
import sys
import glob
import datetime
import time
import subprocess
import job_manager


try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    # DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'
    DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience/butane/dft'


def get_num_species():
    """Function to lookup number of species in the species_list.csv
    """
    species_csv = os.path.join(DFT_DIR, 'species_list.csv')
    species_df = pd.read_csv(species_csv)
    return species_df.i.values[-1]


def index2smiles(species_index):
    """Function to return species smiles given a species index
    looks up the results in the species_list.csv
    """
    species_csv = os.path.join(DFT_DIR, 'species_list.csv')
    species_df = pd.read_csv(species_csv)
    species_index = species_df.SMILES.values[species_index]
    return species_index


def index2name(species_index):
    """Function to return species name given in species_list.csv
    """
    species_csv = os.path.join(DFT_DIR, 'species_list.csv')
    species_df = pd.read_csv(species_csv)
    return species_df[species_df['i'] == species_index]['name'].values[0]


def smiles2index(species_smiles):
    """Function to return species index given a species smiles
    looks up the results in the species_list.csv
    """
    species_csv = os.path.join(DFT_DIR, 'species_list.csv')
    species_df = pd.read_csv(species_csv)
    try:
        species_index = species_df[species_df['SMILES'] == species_smiles]['i'].values[0]
        return species_index
    except IndexError:
        # print(f'could not identify species {species_smiles}')
        raise IndexError(f'could not identify species {species_smiles}')
        # you don't want to equate resonance structures
        import rmgpy.species
        # now we need to check all the species for isomorphism
        ref_sp = rmgpy.species.Species(smiles=species_smiles)
        for i in range(0, len(species_df)):
            sp = rmgpy.species.Species(smiles=species_df['SMILES'].values[i])
            resonance = sp.generate_resonance_structures()
            if resonance:
                sp = resonance
            else:
                sp = [sp]
            for compare_sp in sp:
                if ref_sp.is_isomorphic(compare_sp):
                    return i
        print(f'could not identify species {species_smiles}')


def arkane_complete(species_index):
    """Function to check whether the arkane job is complete for a species
    Expects to find the following directory structure:
    DFT_DIR/thermo/species_XXXX/arkane/RMG_libraries/thermo.py
    Returns True if complete, False otherwise
    """
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    arkane_result = os.path.join(species_dir, 'arkane', 'RMG_libraries', 'thermo.py')
    return os.path.exists(arkane_result)


def termination_status(log_file):
    """Returns:
    0 for Normal termination
    1 for Error termination
    5 for manual skip
    -1 for no termination
    """
    with open(log_file, 'rb') as f:
        f.seek(0, os.SEEK_END)
        normal_termination = False
        error_termination = False
        for i in range(0, 5):
            try:
                f.seek(-2, os.SEEK_CUR)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            saved_position = f.tell()
            last_line = f.readline().decode()
            f.seek(saved_position, os.SEEK_SET)
            if 'Normal termination' in last_line:
                return 0
            elif 'MANUAL SKIP' in last_line.upper():
                return 5
            elif 'Error termination' in last_line:
                return 1
        return -1


def get_n_runs(slurm_array_file):
    """Reads the run.sh file to figure out how many conformers or rotors were meant to run
    """
    with open(slurm_array_file, 'r') as f:
        for line in f:
            if 'SBATCH --array=' in line:
                token = line.split('-')[-1]
                n_runs = 1 + int(token.split('%')[0])
                return n_runs
    return 0


def incomplete_conformers(species_index):
    """Returns a list of indices of incomplete conformers that need to be rerun
    count 'Error termination' as well as 'normal termination'
    Does not work on restart.sh, which has ','
    """
    conformer_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'conformers')

    # Get #conformers from the array job script
    slurm_array_file = os.path.join(conformer_dir, 'run.sh')
    if not os.path.exists(slurm_array_file):
        return True  # no conformers run yet
    n_conformers = get_n_runs(slurm_array_file)

    incomplete_cfs = []
    for cf_index in range(0, n_conformers):
        conformer_file = os.path.join(conformer_dir, f'conformer_{cf_index:04}.log')
        if not os.path.exists(conformer_file):
            incomplete_cfs.append(cf_index)
            continue
        status = termination_status(conformer_file)
        if status == -1:
            incomplete_cfs.append(cf_index)
    return incomplete_cfs


def incomplete_rotors(species_index):
    """Returns a list of indices of incomplete rotors that need to be rerun
    count 'Error termination' as well as 'normal termination'
    Does not work on restart.sh, which has ','
    """
    rotor_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'rotors')

    # Get #rotors from the array job script
    slurm_array_file = os.path.join(rotor_dir, 'run_rotor_calcs.sh')
    if not os.path.exists(slurm_array_file):
        return True  # no rotors run yet
    n_rotors = get_n_runs(slurm_array_file)

    incomplete_rs = []
    for r_index in range(0, n_rotors):
        rotor_file = os.path.join(rotor_dir, f'rotor_{r_index:04}.log')
        if not os.path.exists(rotor_file):
            incomplete_rs.append(r_index)
            continue
        status = termination_status(rotor_file)
        if status == -1:
            incomplete_rs.append(r_index)
    return incomplete_rs


def conformers_complete(species_index):
    """Function to check whether all of the Gaussian conformer jobs have finished running.
    Looks at the run.sh script to find the highest conformer index, then searches each .log file
    for Normal termination
    """
    if incomplete_conformers(species_index):
        return False
    return True


def rotors_complete(species_index):
    """Function to check whether all of the Gaussian rotor jobs have finished running.
    Looks at the run.sh script to find the highest rotor index, then searches each .log file
    for Normal termination
    """
    if incomplete_rotors(species_index):
        return False
    return True


def restart_conformers(species_index):
    """Function to rerun the conformers that didn't converge in time
    """
    # create a new slurm job file to run on west partition, 10 at a time, 2 week max
    missing_conformers = incomplete_conformers(species_index)
    missing_conformers_str = [str(i) for i in missing_conformers]
    indices_str = ','.join(missing_conformers_str)
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    conformer_dir = os.path.join(species_dir, 'conformers')

    # TODO put restart in the gaussian job file
    slurm_run_file = os.path.join(conformer_dir, 'restart.sh')
    slurm_settings = {
        '--job-name': f'g16_cf_{species_index}',
        '--error': 'error.log',
        '--nodes': 1,
        '--partition': 'west',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '14-00:00:00',
        '--cpus-per-task': 16,
        '--array': f'{indices_str}%10',
    }
    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        'fname="conformer_${RUN_i}.com"\n\n',

        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # copy the file and add a restart? this is so messy, but I'm gonna do it
    for cf_idx in missing_conformers:
        pass

        # TODO see if conditions are right to restart in Gaussian:
        # chk file exists
        # previous run made it at least one step in the optimization

    # restart the conformers
    # submit the job
    start_dir = os.getcwd()
    os.chdir(conformer_dir)
    gaussian_conformers_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    gaussian_conformers_job.submit(slurm_cmd)
    os.chdir(start_dir)
    gaussian_conformers_job.wait_all(check_interval=600)


def restart_rotors(species_index):
    """Function to rerun the conformers that didn't converge in time
    """
    # create a new slurm job file to run on west partition, 10 at a time, 2 week max
    missing_rotors = incomplete_rotors(species_index)
    missing_rotors_str = [str(i) for i in missing_rotors]
    indices_str = ','.join(missing_rotors_str)
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    rotor_dir = os.path.join(species_dir, 'rotors')

    # TODO put restart in the gaussian job file
    slurm_run_file = os.path.join(rotor_dir, 'restart.sh')
    slurm_settings = {
        '--job-name': f'g16_rotor_{species_index}',
        '--error': 'error.log',
        '--nodes': 1,
        '--partition': 'west',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '14-00:00:00',
        '--cpus-per-task': 16,
        '--array': f'{indices_str}%10',
    }
    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        'fname="rotor_${RUN_i}.com"\n\n',

        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # submit the job
    start_dir = os.getcwd()
    os.chdir(rotor_dir)
    gaussian_rotors_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    gaussian_rotors_job.submit(slurm_cmd)
    os.chdir(start_dir)
    gaussian_rotors_job.wait_all(check_interval=600)


def run_conformers_job(species_index):
    """Function to call snakemake rule to run conformers
    This function waits until all SLURM jobs are done, so it could take days
    """
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    conformer_dir = os.path.join(species_dir, 'conformers')
    os.makedirs(conformer_dir, exist_ok=True)
    logfile = os.path.join(conformer_dir, 'conformers.log')
    start = time.time()
    timestamp = datetime.datetime.now()
    with open(logfile, 'a') as f:
        f.write(f'Starting conformers job: {timestamp}' + '\n')

    # check if the run was already completed
    if conformers_complete(species_index):
        print('Conformers already ran')
        with open(logfile, 'a') as f:
            f.write('Conformers already ran\n')
        return True

    # TODO make this path relative to the job.py script
    workflow_dir = "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow"

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

    # rerun any conformer jobs that failed to converge in time:
    if not conformers_complete(species_index):
        with open(logfile, 'a') as f:
            f.write('Setting up conformer restart job\n')
        restart_conformers(species_index)  # this waits for jobs to finish
        if not conformers_complete(species_index):
            with open(logfile, 'a') as f:
                f.write('Conformer restart failed\n')
            return False

    end = time.time()
    duration = end - start
    print(f'Gaussian conformer jobs completed in {duration} seconds' + '\n')

    with open(logfile, 'a') as f:
        f.write(f'Gaussian conformer jobs completed in {duration} seconds' + '\n')

    return True


def read_gaussian_energy(logfile):
    with open(logfile, 'r') as f:
        for line in f:
            if 'Sum of electronic and zero-point Energies= ' in line:
                energy = float(line.split()[-1])
                return energy
    return 0


def get_lowest_conformer(species_index):
    """Returns the filepath of the lowest energy conformer logfile
    """
    conformer_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}', 'conformers')
    slurm_array_file = os.path.join(conformer_dir, 'run.sh')
    if not os.path.exists(slurm_array_file):
        return None  # no conformers run yet
    n_conformers = get_n_runs(slurm_array_file)
    lowest_energy = 999999
    best_conformer_file = None
    for cf_index in range(0, n_conformers):
        conformer_file = os.path.join(conformer_dir, f'conformer_{cf_index:04}.log')
        status = termination_status(conformer_file)
        if status != 0:
            continue
        energy = read_gaussian_energy(conformer_file)
        print(cf_index, energy)
        if energy < lowest_energy:
            lowest_energy = energy
            best_conformer_file = conformer_file

    return best_conformer_file


def run_rotors_job(species_index):
    # start a job that calls snakemake to run rotors

    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    rotor_dir = os.path.join(species_dir, 'rotors')
    os.makedirs(rotor_dir, exist_ok=True)
    logfile = os.path.join(rotor_dir, 'rotors.log')
    start = time.time()
    timestamp = datetime.datetime.now()
    with open(logfile, 'a') as f:
        f.write(f'Starting rotors job: {timestamp}' + '\n')

    # check if a rotor job was already completed
    if os.path.exists(os.path.join(rotor_dir, 'NO_ROTORS.txt')):
        print('No rotors to run')
        with open(logfile, 'a') as f:
            f.write('No rotors to run\n')
        return True
    elif rotors_complete(species_index):
        print('Rotors already ran')
        with open(logfile, 'a') as f:
            f.write('Rotors already ran\n')
        return True

    rotor_cmd = f'snakemake -c1 run_rotors --config species_index={species_index}'
    print(f'Running {rotor_cmd}')
    cmd_pieces = rotor_cmd.split()
    proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
    print(proc)

    # wait 5 minutes for the rotor gaussian job to begin
    time.sleep(300)
    g16_job_number = ''

    rotor_slurm_files = glob.glob(os.path.join(rotor_dir, 'slurm-*'))
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
    rotor_job.wait_all(check_interval=30)

    # rerun any rotor jobs that failed to converge in time:
    if not rotors_complete(species_index):
        with open(logfile, 'a') as f:
            f.write('Setting up rotor restart job\n')
        restart_rotors(species_index)  # this waits for jobs to finish
        if not rotors_complete(species_index):
            with open(logfile, 'a') as f:
                f.write('Rotor restart failed\n')
            return False

    end = time.time()
    duration = end - start
    print(f'Gaussian rotor jobs completed in {duration} seconds' + '\n')
    with open(logfile, 'a') as f:
        f.write(f'Gaussian rotor jobs completed in {duration} seconds' + '\n')
    return True


def run_arkane_job(species_index):
    # start a job that calls snakemake to run arkane
    species_dir = os.path.join(DFT_DIR, 'thermo', f'species_{species_index:04}')
    arkane_dir = os.path.join(species_dir, 'arkane')
    os.makedirs(arkane_dir, exist_ok=True)
    arkane_result = os.path.join(arkane_dir, 'RMG_libraries', 'thermo.py')
    if arkane_complete(species_index):
        print('Arkane job already ran')
        return True

    cwd = os.getcwd()
    snakemake_dir = '/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow'
    os.chdir(snakemake_dir)
    arkane_cmd = f'snakemake -c1 run_arkane_thermo --config species_index={species_index}'
    print(f'Running {arkane_cmd}')
    cmd_pieces = arkane_cmd.split()
    proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
    print(proc)
    os.chdir(cwd)

    # wait 10 minutes for Arkane start/finish
    # try to read the slurm file in
    print('Waiting for arkane job')
    logfile = os.path.join(arkane_dir, 'snakemake_arkane.log')
    with open(logfile, 'a') as f:
        f.write('Waiting for arkane job\n')
    while not os.path.exists(arkane_result):
        time.sleep(30)

        # TODO, give up if it has started running but hasn't completed in twenty minutes
    print('Arkane complete')
    with open(logfile, 'a') as f:
        f.write('Arkane complete\n')

    end = time.time()
    duration = end - start
    print(f'COMPLETED {species_smiles} IN {duration} SECONDS')
    with open(logfile, 'a') as f:
        f.write(f'COMPLETED {species_smiles} IN {duration} SECONDS' + '\n')


# temporary function to make no_rotors library -- delete this after you're done with it
def run_arkane_job_no_rotors(species_index):
    # start a job that calls snakemake to run arkane
    species_dir = os.path.join(DFT_DIR, 'no_rotors_thermo', f'species_{species_index:04}')
    arkane_dir = os.path.join(species_dir, 'arkane')
    os.makedirs(arkane_dir, exist_ok=True)
    arkane_result = os.path.join(arkane_dir, 'RMG_libraries', 'thermo.py')
    if os.path.exists(f'/work/westgroup/harris.se/autoscience/autoscience/butane/dft/no_rotors_thermo/species_{species_index:04}/arkane/RMG_libraries/thermo.py'):
        print("arkane already ran")
        return

    cwd = os.getcwd()
    snakemake_dir = '/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow'
    os.chdir(snakemake_dir)
    arkane_cmd = f'snakemake -c1 run_arkane_thermo --config species_index={species_index}'
    print(f'Running {arkane_cmd}')
    cmd_pieces = arkane_cmd.split()
    proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
    print(proc)
    os.chdir(cwd)

    # wait 10 minutes for Arkane start/finish
    # try to read the slurm file in
    print('Waiting for arkane job')
    logfile = os.path.join(arkane_dir, 'snakemake_arkane.log')
    with open(logfile, 'a') as f:
        f.write('Waiting for arkane job\n')
    while not os.path.exists(arkane_result):
        time.sleep(30)

        # TODO, give up if it has started running but hasn't completed in twenty minutes
    print('Arkane complete')
    with open(logfile, 'a') as f:
        f.write('Arkane complete\n')

    end = time.time()
    duration = end - start
    print(f'COMPLETED {species_smiles} IN {duration} SECONDS')
    with open(logfile, 'a') as f:
        f.write(f'COMPLETED {species_smiles} IN {duration} SECONDS' + '\n')
