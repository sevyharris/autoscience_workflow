# Functions for running a kinetics job using this workflow
import pandas as pd
import numpy as np
import re
import os
import glob
import numpy as np
import subprocess
import autotst.reaction
import autotst.calculator.vibrational_analysis
import autotst.calculator.gaussian
try:
    import job_manager
except ImportError:
    pass
try:
    from hotbit import Hotbit  # TODO - move this and other autoTST dependencies elsewhere
except ImportError:
    pass
import ase.io.gaussian
import ase.calculators.lj
import rmgpy.reaction
import rmgpy.species




try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    # DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'
    DFT_DIR = "/work/westgroup/harris.se/autoscience/autoscience/butane/dft"


def ordered_array_str(list_of_indices):
    # convenient script for putting a list of task numbers into a string that can be used for a SLURM array job
    # assume it's sorted
    if len(list_of_indices) == 1:
        return str(list_of_indices[0])
    elif len(list_of_indices) == 2:
        return f'{list_of_indices[0]},{list_of_indices[1]}'
    elif len(list_of_indices) == 0:
        return -1

    array_str = str(list_of_indices[0]) + '-'
    for j in range(1, len(list_of_indices) - 1):
        if list_of_indices[j] - list_of_indices[j - 1] != 1:
            if j > 1:
                array_str += str(list_of_indices[j - 1])
                array_str += f',{list_of_indices[j]}-'
            else:
                array_str = array_str.replace('-', ',')
                array_str += f'{list_of_indices[j]}-'
        # cap the end
        if j + 2 == len(list_of_indices):
            if list_of_indices[j + 1] - list_of_indices[j] != 1:
                array_str += f'{list_of_indices[j]},{list_of_indices[j + 1]}'
            else:
                array_str += f'{list_of_indices[j + 1]}'

    return array_str


def get_num_reactions():
    """Function to lookup number of reactions in the reaction_list.csv
    """
    reaction_csv = os.path.join(DFT_DIR, 'reaction_list.csv')
    reaction_df = pd.read_csv(reaction_csv)
    return reaction_df.i.values[-1]


def reaction2smiles(reaction):
    """Takes an RMG reaction and returns the smiles representation
    """
    string = ""
    for react in reaction.reactants:
        if isinstance(react, rmgpy.species.Species):
            string += f"{react.molecule[0].to_smiles()}+"
        elif isinstance(react, rmgpy.molecule.Molecule):
            string += f"{react.to_smiles()}+"
    string = string[:-1]
    string += "_"
    for prod in reaction.products:
        if isinstance(prod, rmgpy.species.Species):
            string += f"{prod.molecule[0].to_smiles()}+"
        elif isinstance(prod, rmgpy.molecule.Molecule):
            string += f"{prod.to_smiles()}+"
    label = string[:-1]
    return label


def smiles2reaction(reaction_smiles):
    """Takes the reaction smiles and produces a corresponding rmg reaction
    """
    reaction = rmgpy.reaction.Reaction()
    reactants = []
    products = []

    # handle CO case
    if '[C-]#[O+]' in reaction_smiles:
        CO = rmgpy.species.Species(smiles='[C-]#[O+]')
        reaction_smiles = reaction_smiles.replace('[C-]#[O+]', 'carbonmonoxide')
    if '[O-][N+]#C' in reaction_smiles:
        CHNO = rmgpy.species.Species(smiles='[O-][N+]#C')
        reaction_smiles = reaction_smiles.replace('[O-][N+]#C', 'formonitrileoxide')
    if '[O-][N+]=C' in reaction_smiles:
        CH2NO = rmgpy.species.Species(smiles='[O-][N+]=C')
        reaction_smiles = reaction_smiles.replace('[O-][N+]=C', 'methylenenitroxide')

    reactant_token = reaction_smiles.split('_')[0]
    product_token = reaction_smiles.split('_')[1]

    reactant_tokens = reactant_token.split('+')
    product_tokens = product_token.split('+')

    # print(product_tokens)
    for reactant_str in reactant_tokens:
        if reactant_str == 'carbonmonoxide':
            reactant_str = '[C-]#[O+]'
        elif reactant_str == 'formonitrileoxide':
            reactant_str = '[O-][N+]#C'
        elif reactant_str == 'methylenenitroxide':
            reactant_str = '[O-][N+]=C'
        reactant = rmgpy.species.Species(smiles=reactant_str)
        reactants.append(reactant)
    
    for product_str in product_tokens:
        if product_str == 'carbonmonoxide':
            product_str = '[C-]#[O+]'
        elif product_str == 'formonitrileoxide':
            product_str = '[O-][N+]#C'
        elif product_str == 'methylenenitroxide':
            ptoduct_str = '[O-][N+]=C'
        # print(product_str)
        product = rmgpy.species.Species(smiles=product_str)
        products.append(product)

    reaction.reactants = reactants
    reaction.products = products
    return reaction


def reaction_index2smiles(reaction_index):
    """Function to return reaction smiles given a reaction index
    looks up the results in the reaction_list.csv
    """
    reaction_csv = os.path.join(DFT_DIR, 'reaction_list.csv')
    reaction_df = pd.read_csv(reaction_csv)
    reaction_smiles = reaction_df['SMILES'].values[reaction_index]
    return reaction_smiles


def reaction_smiles2index(reaction_smiles):
    """Function to return reaction index given a smiles reaction
    doesn't necessarily have to be in the right order
    RMG reaction will check for isomorphism
    """
    # first check to see if the exact smiles is in the CSV
    reaction_csv = os.path.join(DFT_DIR, 'reaction_list.csv')
    reaction_df = pd.read_csv(reaction_csv)
    if reaction_smiles in reaction_df['SMILES'].values:
        idx = np.where(reaction_df['SMILES'].values == reaction_smiles)[0][0]
        return reaction_df['i'].values[idx]
    else:
        # use rmgpy.reaction to check for isomorphism
        ref_reaction = smiles2reaction(reaction_smiles)
        for i in range(0, len(reaction_df)):
            csv_reaction = smiles2reaction(reaction_df['SMILES'].values[i])
            if ref_reaction.is_isomorphic(csv_reaction):
                return i
    # reaction not found
    return -1


def termination_status(log_file):
    """Returns:
    0 for Normal termination
    1 for Error termination not covered below
    2 for Error termination - due to all degrees of freedom being frozen
    3 for Error termination - Problem with the distance matrix.
    4 for No NMR shielding tensors so no spin-rotation constants  # TODO debug this instead of ignoring it
    -1 for no termination
    """
    error_termination = False
    with open(log_file, 'rb') as f:
        f.seek(0, os.SEEK_END)
        normal_termination = False
        error_termination = False
        for i in range(0, 20):
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
            elif 'Error termination' in last_line:
                error_termination = True
            elif 'All variables have been frozen' in last_line:
                return 2
            elif 'Problem with the distance matrix' in last_line:
                return 3
            elif 'No NMR shielding tensors so no spin-rotation constants' in last_line:
                return 4
            elif 'MANUAL SKIP' in last_line.upper():
                return 5
        if error_termination:
            return 1
        return -1


def shell_complete(reaction_index, use_reverse=False):
    # return True if all shell calculations have completed
    # returns true even if one of the shell calculations had the error with distance matrix problem
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{int(reaction_index):04}')
    shell_dir = os.path.join(reaction_base_dir, 'shell')
    if not os.path.exists(shell_dir):
        return False
    logfile = os.path.join(shell_dir, 'shell.log')

    # Check for already finished shell logfiles
    # first, return if all of them have finished
    shell_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        shell_label = 'rev_ts_0000.log'
        direction = 'reverse'
    shell_opt_log = os.path.join(shell_dir, shell_label)

    shell_gaussian_logs = glob.glob(os.path.join(shell_dir, shell_label[:-8] + '*.log'))
    incomplete_indices = []
    good_runs = []
    for shell_opt_log in shell_gaussian_logs:
        if os.path.exists(shell_opt_log):
            status = termination_status(shell_opt_log)
            print(status, shell_opt_log)
            if status == 0:
                matches = re.search(shell_label[:-8] + '([0-9]{4})', shell_opt_log)
                run_index = int(matches[1])
                good_runs.append(run_index)
            elif status == 2 or status == 3 or status == 4 or status == 5:
                # print('Shell optimization already ran')
                with open(logfile, 'a') as f:
                    f.write('Shell optimization already ran\n')
            else:
                matches = re.search(shell_label[:-8] + '([0-9]{4})', shell_opt_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)
    if not incomplete_indices and len(good_runs) > 0:
        return True  # only if there's at least one useable run
    return False


def run_TS_shell_calc(reaction_index, use_reverse=False, max_combos=300, max_conformers=12):
    """Start an optimization keeping the reaction center fixed
    """
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    os.makedirs(reaction_base_dir, exist_ok=True)
    shell_dir = os.path.join(reaction_base_dir, 'shell')
    os.makedirs(shell_dir, exist_ok=True)
    logfile = os.path.join(shell_dir, 'shell.log')

    # Check for already finished shell logfiles
    # first, return if all of them have finished
    shell_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        shell_label = 'rev_ts_0000.log'
        direction = 'reverse'
    shell_opt_log = os.path.join(shell_dir, shell_label)

    shell_gaussian_logs = glob.glob(os.path.join(shell_dir, shell_label[:-8] + '*.log'))
    incomplete_indices = []
    for shell_opt_log in shell_gaussian_logs:
        if os.path.exists(shell_opt_log):
            status = termination_status(shell_opt_log)
            if status == 0 or status == 2 or status == 3 or status == 4 or status == 5:
                print('Shell optimization already ran')
                with open(logfile, 'a') as f:
                    f.write('Shell optimization already ran\n')
            else:
                matches = re.search(shell_label[:-8] + '([0-9]{4})', shell_opt_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)
    if not incomplete_indices and len(shell_gaussian_logs) > 0:
        return True  # only if all of them ran

    print('Constructing reaction in AutoTST...')
    with open(logfile, 'a') as f:
        f.write('Constructing reaction in AutoTST...\n')
    reaction_smiles = reaction_index2smiles(reaction_index)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    reaction.ts[direction][0].get_molecules()
    reaction.generate_conformers(ase_calculator=Hotbit(), max_combos=max_combos, max_conformers=max_conformers)
    print('Done generating conformers in AutoTST...')
    print(f'{len(reaction.ts[direction])} conformers found')
    with open(logfile, 'a') as f:
        f.write('Done generating conformers in AutoTST...\n')
        f.write(f'{len(reaction.ts[direction])} conformers found' + '\n')

    # Do the shell calculations
    # write Gaussian input files
    slurm_array_idx = []
    restart = False
    for i in range(0, len(reaction.ts[direction])):
        if i not in incomplete_indices and len(shell_gaussian_logs) > 0:
            with open(logfile, 'a') as f:
                f.write(f'skipping completed shell {i}' + '\n')
            restart = True
            continue

        shell_label = shell_label[:-8] + f'{i:04}.log'
        # TODO if the thing already ran, do not rerun it
        slurm_array_idx.append(i)

        ts = reaction.ts[direction][i]
        gaussian = autotst.calculator.gaussian.Gaussian(conformer=ts)
        calc = gaussian.get_shell_calc()
        calc.label = shell_label[:-4]
        calc.directory = shell_dir
        calc.parameters.pop('scratch')
        calc.parameters.pop('multiplicity')
        calc.parameters['mult'] = ts.rmg_molecule.multiplicity
        calc.write_input(ts.ase_molecule)

        # Get rid of double-space between xyz block and mod-redundant section
        double_space = False
        lines = []
        with open(os.path.join(shell_dir, calc.label + '.com'), 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                if lines[j] == '\n' and lines[j - 1] == '\n':
                    double_space = True
                    break
        if double_space:
            lines = lines[0:j - 1] + lines[j:]
            with open(os.path.join(shell_dir, calc.label + '.com'), 'w') as f:
                f.writelines(lines)

    # make the shell slurm script
    slurm_run_file = os.path.join(shell_dir, f'run_shell_opt.sh')
    slurm_settings = {
        '--job-name': f'g16_shell_{reaction_index}',
        '--error': 'error.log',
        '--output': 'output.log',
        '--nodes': 1,
        '--partition': 'west,short',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '24:00:00',
        '--cpus-per-task': 16,
        '--array': ordered_array_str(slurm_array_idx),
        # TODO make this an array if multiple forward ts's
    }
    if restart:
        slurm_run_file = os.path.join(shell_dir, f'restart.sh')
        slurm_settings = {
            '--job-name': f'g16_shell_restart_{reaction_index}',
            '--error': 'restart_error.log',
            '--output': 'restart_output.log',
            '--nodes': 1,
            '--partition': 'short',
            '--mem': '20Gb',
            '--time': '24:00:00',
            '--cpus-per-task': 32,
            '--array': ordered_array_str(slurm_array_idx),
        }

    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        f'fname="{shell_label[:-8]}' + '${RUN_i}.com"\n\n',

        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # submit the job
    with open(logfile, 'a') as f:
        f.write('Submitting shell optimization job\n')
    start_dir = os.getcwd()
    os.chdir(shell_dir)
    shell_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    shell_job.submit(slurm_cmd)

    # only wait after all jobs have been submitted
    os.chdir(start_dir)
    shell_job.wait(check_interval=60)


def run_TS_center_calc(reaction_index, use_reverse=False):
    """Start a constrained saddle search with the reaction center fixed
    """
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    shell_dir = os.path.join(reaction_base_dir, 'shell')
    center_dir = os.path.join(reaction_base_dir, 'center')
    os.makedirs(center_dir, exist_ok=True)
    logfile = os.path.join(center_dir, 'center.log')

    # Check for already finished center logfiles
    # first, return if all of them have finished
    shell_label = 'fwd_ts_0000.log'
    center_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        shell_label = 'rev_ts_0000.log'
        center_label = 'rev_ts_0000.log'
        direction = 'reverse'
    shell_opt_log = os.path.join(shell_dir, shell_label)
    center_opt_log = os.path.join(center_dir, center_label)

    shell_gaussian_logs = glob.glob(os.path.join(shell_dir, shell_label[:-8] + '*.log'))
    center_gaussian_logs = glob.glob(os.path.join(center_dir, center_label[:-8] + '*.log'))
    incomplete_indices = []
    for center_opt_log in center_gaussian_logs:
        if os.path.exists(center_opt_log):
            status = termination_status(center_opt_log)
            if status == 0 or status == 2 or status == 3 or status == 4 or status == 5:
                print('Center saddle point already ran')
                with open(logfile, 'a') as f:
                    f.write('Center saddle point already ran\n')
            else:
                matches = re.search(shell_label[:-8] + '([0-9]{4})', center_opt_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)
    if not incomplete_indices and len(center_gaussian_logs) > 0:
        return True  # only if all of them ran

    print('Constructing reaction in AutoTST...')
    with open(logfile, 'a') as f:
        f.write('Constructing reaction in AutoTST...\n')
    reaction_smiles = reaction_index2smiles(reaction_index)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    reaction.ts[direction][0].get_molecules()
    reaction.generate_conformers(ase_calculator=Hotbit())
    print('Done generating conformers in AutoTST...')
    print(f'{len(reaction.ts[direction])} conformers found')
    with open(logfile, 'a') as f:
        f.write('Done generating conformers in AutoTST...\n')
        f.write(f'{len(reaction.ts[direction])} conformers found' + '\n')
    # TODO - a way to export the reaction conformers from the shell run so we don't have to repeat it here?

    # define incomplete indices
    center_gaussian_logs = glob.glob(os.path.join(center_dir, center_label[:-8] + '*.log'))
    incomplete_indices = []
    for center_ts_log in center_gaussian_logs:
        if os.path.exists(center_ts_log):
            status = termination_status(center_ts_log)
            if status == -1 or status == 1:
                matches = re.search(center_label[:-8] + '([0-9]{4})', center_ts_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)

    restart = False
    slurm_array_idx = []
    for i in range(0, len(reaction.ts[direction])):
        if i not in incomplete_indices and len(center_gaussian_logs) > 0:
            print(f'Skipping completed index {i}')
            restart = True
            continue

        center_label = center_label[:-8] + f'{i:04}.log'
        shell_opt = os.path.join(shell_dir, center_label)

        # skip shell conformers that didn't converge
        status = termination_status(shell_opt)
        if status != 0:
            print(f'Skipping unconverged shell {i}')
            with open(logfile, 'a') as f:
                f.write(f'Skipping unconverged shell {i}')
            continue

        try:
            with open(shell_opt, 'r') as f:
                atoms = ase.io.gaussian.read_gaussian_out(f)
                reaction.ts[direction][i]._ase_molecule = atoms
        except IndexError:
            # handle case where all degrees of freedom were frozen in the shell calculation
            if len(reaction.ts[direction][i]._ase_molecule) > 3:
                raise ValueError('Shell optimization failed to converge. Rerun it!')

        slurm_array_idx.append(i)

        ts = reaction.ts[direction][i]
        gaussian = autotst.calculator.gaussian.Gaussian(conformer=ts)
        # calc = gaussian.get_overall_calc()
        calc = gaussian.get_center_calc()
        calc.label = center_label[:-4]
        calc.directory = center_dir
        calc.parameters.pop('scratch')
        calc.parameters.pop('multiplicity')
        calc.parameters['mult'] = ts.rmg_molecule.multiplicity
        calc.write_input(ts.ase_molecule)

        # Get rid of double-space between xyz block and mod-redundant section
        double_space = False
        lines = []
        with open(os.path.join(center_dir, calc.label + '.com'), 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                if lines[j] == '\n' and lines[j - 1] == '\n':
                    double_space = True
                    break
        if double_space:
            lines = lines[0:j - 1] + lines[j:]
            with open(os.path.join(center_dir, calc.label + '.com'), 'w') as f:
                f.writelines(lines)

    # make the overall slurm script
    slurm_run_file = os.path.join(center_dir, f'run_center_opt.sh')
    if restart:
        slurm_run_file = os.path.join(center_dir, f'restart.sh')
    slurm_settings = {
        '--job-name': f'g16_center_{reaction_index}',
        '--error': 'error.log',
        '--output': 'output.log',
        '--nodes': 1,
        '--partition': 'west,short',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '24:00:00',
        '--cpus-per-task': 16,
        '--array': ordered_array_str(slurm_array_idx),
    }

    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        f'fname="{center_label[:-8]}' + '${RUN_i}.com"\n\n',

        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # submit the job
    start_dir = os.getcwd()
    os.chdir(center_dir)
    center_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    center_job.submit(slurm_cmd)

    # only wait once all jobs have been submitted
    os.chdir(start_dir)
    center_job.wait(check_interval=60)


def overall_complete(reaction_index, use_reverse=False):
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{int(reaction_index):04}')
    os.makedirs(reaction_base_dir, exist_ok=True)
    overall_dir = os.path.join(reaction_base_dir, 'overall')

    logfile = os.path.join(overall_dir, 'overall.log')

    overall_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        overall_label = 'rev_ts_0000.log'
        direction = 'reverse'
    overall_ts_log = os.path.join(overall_dir, overall_label)

    # check that the overall job hasn't finished
    # in this case, terminating with an error is okay, as long as it terminated
    overall_gaussian_logs = glob.glob(os.path.join(overall_dir, overall_label[:-8] + '*.log'))
    incomplete_indices = []
    for overall_ts_log in overall_gaussian_logs:
        if os.path.exists(overall_ts_log):
            status = termination_status(overall_ts_log)
            print(status, overall_ts_log)
            if status == 0:
                print(f'Overall TS optimization already ran for reaction {overall_ts_log}')
                with open(logfile, 'a') as f:
                    f.write(f'Overall TS optimization already ran for reaction {overall_ts_log}' + '\n')
            elif status == 2 or status == 3 or status == 4 or status == 5:  # completed but with error
                matches = re.search(overall_label[:-8] + '([0-9]{4})', overall_ts_log)
                run_index = int(matches[1])
                # incomplete_indices.append(run_index)
                print(f'Overall TS optimization ran with error for {overall_ts_log}')
                with open(logfile, 'a') as f:
                    f.write(f'Overall TS optimization ran with error for {overall_ts_log}' + '\n')
            else:
                matches = re.search(overall_label[:-8] + '([0-9]{4})', overall_ts_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)
    print(incomplete_indices)
    if not incomplete_indices and len(overall_gaussian_logs) > 0:
        return True
    return False


def run_TS_overall_calc(reaction_index, use_reverse=False, max_combos=300, max_conformers=12):
    """Start a TS optimization from the geometry of the shell calculation
    """
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    os.makedirs(reaction_base_dir, exist_ok=True)

    shell_dir = os.path.join(reaction_base_dir, 'shell')
    overall_dir = os.path.join(reaction_base_dir, 'overall')
    os.makedirs(overall_dir, exist_ok=True)

    logfile = os.path.join(overall_dir, 'overall.log')

    overall_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        overall_label = 'rev_ts_0000.log'
        direction = 'reverse'
    overall_ts_log = os.path.join(overall_dir, overall_label)

    # check that the overall job hasn't finished
    # in this case, terminating with an error is okay, as long as it terminated (status 2)
    if overall_complete(reaction_index):
        return True

    print('Constructing reaction in AutoTST...')
    with open(logfile, 'a') as f:
        f.write('Constructing reaction in AutoTST...\n')
    reaction_smiles = reaction_index2smiles(reaction_index)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    reaction.ts[direction][0].get_molecules()
    reaction.generate_conformers(ase_calculator=Hotbit(), max_combos=max_combos, max_conformers=max_conformers)
    print('Done generating conformers in AutoTST...')
    print(f'{len(reaction.ts[direction])} conformers found')
    with open(logfile, 'a') as f:
        f.write('Done generating conformers in AutoTST...\n')
        f.write(f'{len(reaction.ts[direction])} conformers found' + '\n')
    # TODO - a way to export the reaction conformers from the shell run so we don't have to repeat it here?

    # define incomplete indices
    overall_gaussian_logs = glob.glob(os.path.join(overall_dir, overall_label[:-8] + '*.log'))
    incomplete_indices = []
    for overall_ts_log in overall_gaussian_logs:
        if os.path.exists(overall_ts_log):
            status = termination_status(overall_ts_log)
            if status == -1 or status == 1:
                matches = re.search(overall_label[:-8] + '([0-9]{4})', overall_ts_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)

    restart = False
    slurm_array_idx = []
    for i in range(0, len(reaction.ts[direction])):
        if i not in incomplete_indices and len(overall_gaussian_logs) > 0:
            print(f'Skipping completed index {i}')
            restart = True
            continue

        overall_label = overall_label[:-8] + f'{i:04}.log'
        shell_opt = os.path.join(shell_dir, overall_label)
        
        if not os.path.exists(shell_opt):
            print(f'WHY does it think this shell log should exist??? {shell_opt}')
            with open(logfile, 'a') as f:
                f.write(f'WHY does it think this shell log should exist??? {shell_opt}')
            continue

        # skip shell conformers that didn't converge
        status = termination_status(shell_opt)
        if status != 0:
            print(f'Skipping unconverged shell {i}')
            with open(logfile, 'a') as f:
                f.write(f'Skipping unconverged shell {i}')
            continue

        try:
            with open(shell_opt, 'r') as f:
                atoms = ase.io.gaussian.read_gaussian_out(f)
                reaction.ts[direction][i]._ase_molecule = atoms
        except IndexError:
            # handle case where all degrees of freedom were frozen in the shell calculation
            if len(reaction.ts[direction][i]._ase_molecule) > 3:
                raise ValueError('Shell optimization failed to converge. Rerun it!')

        slurm_array_idx.append(i)

        ts = reaction.ts[direction][i]
        gaussian = autotst.calculator.gaussian.Gaussian(conformer=ts)
        calc = gaussian.get_overall_calc()
        calc.label = overall_label[:-4]
        calc.directory = overall_dir
        calc.parameters.pop('scratch')
        calc.parameters.pop('multiplicity')
        calc.parameters['mult'] = ts.rmg_molecule.multiplicity
        calc.write_input(ts.ase_molecule)

        # Get rid of double-space between xyz block and mod-redundant section
        double_space = False
        lines = []
        with open(os.path.join(overall_dir, calc.label + '.com'), 'r') as f:
            lines = f.readlines()
            for j in range(1, len(lines)):
                if lines[j] == '\n' and lines[j - 1] == '\n':
                    double_space = True
                    break
        if double_space:
            lines = lines[0:j - 1] + lines[j:]
            with open(os.path.join(overall_dir, calc.label + '.com'), 'w') as f:
                f.writelines(lines)

    # make the overall slurm script
    slurm_run_file = os.path.join(overall_dir, f'run_overall_opt.sh')
    if restart:
        slurm_run_file = os.path.join(overall_dir, f'restart.sh')
    slurm_settings = {
        '--job-name': f'g16_overall_{reaction_index}',
        '--error': 'error.log',
        '--output': 'output.log',
        '--nodes': 1,
        '--partition': 'west,short',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '24:00:00',
        '--cpus-per-task': 16,
        '--array': ordered_array_str(slurm_array_idx),
    }

    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        f'fname="{overall_label[:-8]}' + '${RUN_i}.com"\n\n',

        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # submit the job
    start_dir = os.getcwd()
    os.chdir(overall_dir)
    overall_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    overall_job.submit(slurm_cmd)

    # only wait once all jobs have been submitted
    os.chdir(start_dir)
    overall_job.wait(check_interval=60)


def arkane_complete(reaction_index):
    return os.path.exists(os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}', 'arkane', 'RMG_libraries', 'reactions.py'))


def run_arkane_job(reaction_index):
    # start a job to run arkane
    reaction_smiles = reaction_index2smiles(reaction_index)
    print(f'starting run_arkane_job for reaction {reaction_index} {reaction_smiles}')
    species_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    arkane_dir = os.path.join(species_dir, 'arkane')
    os.makedirs(arkane_dir, exist_ok=True)
    arkane_result = os.path.join(arkane_dir, 'RMG_libraries', 'reactions.py')  # ??
    if arkane_complete(reaction_index):
        print('Arkane job already ran')
        return True

    # setup the arkane job
    # arkane_cmd = f'snakemake -c1 run_arkane_kinetics --config reaction_index={reaction_index}'
    # print(f'Running {arkane_cmd}')
    # cmd_pieces = arkane_cmd.split()
    # proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=None, stderr=None, close_fds=True)
    # print(proc)

    # setup the arkane job using the script directly since conda environments are incompatible
    setup_script = "setup_arkane_kinetics.py"
    setup_script = "/work/westgroup/harris.se/autoscience/autoscience_workflow/workflow/scripts/kinetics/setup_arkane_kinetics.py"
    setup_cmd = f'python {setup_script} {reaction_index}'
    cmd_pieces = setup_cmd.split()
    # apparently subprocess.call is blocking and subprocess.Popen is not
    proc = subprocess.call(cmd_pieces)

    # Run the arkane job
    start_dir = os.getcwd()
    os.chdir(arkane_dir)
    # arkane_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch run_arkane.sh"
    slurm_pieces = slurm_cmd.split()
    proc = subprocess.call(slurm_pieces)
    # arkane_job.submit(slurm_cmd)


def run_vibrational_analysis(reaction_smiles, reaction_logfile):
    # runs the vibrational analysis check for a given TS, returns True if the TS is confirmed
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    va = autotst.calculator.vibrational_analysis.VibrationalAnalysis(
        transitionstate=reaction.ts['forward'][0], log_file=reaction_logfile
    )
    result = va.validate_ts()
    return result


def vibrational_analysis_confirms_ts(reaction_index):
    # Check whether TS is confirmed by vibrational analysis alone
    # If confirmed, there will be a vibrational_analysis_check.txt file in the arkane directory with True
    reaction_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    vib_file = os.path.join(reaction_dir, 'arkane', 'vibrational_analysis_check.txt')

    if not os.path.exists(vib_file):
        return False

    with open(vib_file, 'r') as f:
        vib_check = f.read()
    if vib_check == 'True':
        return True
    return False


def run_IRC_check(reaction_index):
    # TODO get this to run using only smiles
    reaction_smiles = reaction_index2smiles(reaction_index)
    print(f'starting run_IRC_check for reaction {reaction_index} {reaction_smiles}')
    reaction_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    irc_dir = os.path.join(reaction_dir, 'irc')
    logfile = os.path.join(irc_dir, 'irc.log')
    os.makedirs(irc_dir, exist_ok=True)

    # Try the vibrational analysis check first
    reaction_logfiles = glob.glob(
        os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}', 'arkane', 'fwd_*.log')
    )
    assert len(reaction_logfiles) == 1
    reaction_logfile = reaction_logfiles[0]

    vib_check_result = run_vibrational_analysis(reaction_smiles, reaction_logfile)
    # save result to file
    # with open(os.path.join(irc_dir, 'vibrational_analysis_check.txt'), 'w') as f:  # maybe it does belong in irc folder
    with open(os.path.join(reaction_dir, 'arkane', 'vibrational_analysis_check.txt'), 'w') as f:
        f.write(str(vib_check_result))

    if vib_check_result:
        # IRC check is confirmed with vibrational analysis
        return vib_check_result

    # If vibrational analysis check fails, run full IRC check using geometry from reaction_logfile
    # setup the IRC job
    # TODO check for previous run of IRC

    # figure out if we need to restart the IRC
    if os.path.exists(os.path.join(irc_dir, 'irc_result.txt')):
        with open(irc_dir, 'irc_result.txt', 'r') as f:
            irc_result = f.read()
        if irc_result == 'True':
            return True
        else:
            # TODO - restart the IRC
            pass

    print('Constructing reaction in AutoTST...')
    with open(logfile, 'a') as f:
        f.write('Constructing reaction in AutoTST...\n')

    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    direction = 'forward'
    reaction.ts[direction][0].get_molecules()
    # reaction.generate_conformers(ase_calculator=Hotbit())
    # reaction.generate_conformers(ase_calculator=ase.calculators.lj.LennardJones())
    reaction.generate_conformers(ase_calculator='SKIP')

    print('Done generating conformers in AutoTST...')
    print(f'{len(reaction.ts[direction])} conformers found')
    with open(logfile, 'a') as f:
        f.write('Done generating conformers in AutoTST...\n')
        f.write(f'{len(reaction.ts[direction])} conformers found' + '\n')
    # TODO - a way to export the reaction conformers from the shell run so we don't have to repeat it here?

    # get the conformer index associated with the reaction log file
    conformer_index = int(reaction_logfile.split('_')[-1].split('.')[0])
    assert conformer_index < len(reaction.ts[direction])
    irc_label = f'fwd_ts_{conformer_index:04}.log'

    # read in geometry from reaction log file
    with open(reaction_logfile, 'r') as f:
        atoms = ase.io.gaussian.read_gaussian_out(f)
        reaction.ts[direction][conformer_index]._ase_molecule = atoms

    ts = reaction.ts[direction][conformer_index]
    gaussian = autotst.calculator.gaussian.Gaussian(conformer=ts)
    calc = gaussian.get_irc_calc()
    calc.label = irc_label[:-4]
    calc.directory = irc_dir
    calc.parameters.pop('scratch')
    calc.parameters.pop('multiplicity')
    calc.parameters['mult'] = ts.rmg_molecule.multiplicity
    calc.write_input(ts.ase_molecule)

    # make the overall slurm script
    slurm_run_file = os.path.join(irc_dir, f'run_irc.sh')
    slurm_settings = {
        '--job-name': f'g16_irc_{reaction_index}',
        '--error': 'error.log',
        '--output': 'output.log',
        '--nodes': 1,
        '--partition': 'west,short',
        '--exclude': 'c5003',
        '--mem': '20Gb',
        '--time': '24:00:00',
        '--cpus-per-task': 16,
    }

    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/gaussian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        # 'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        # f'fname="{overall_label[:-8]}' + '${RUN_i}.com"\n\n',

        f'fname="{irc_label[:-4]}.com"\n\n',
        'g16 $fname\n',
    ]
    slurm_file_writer.write_file()

    # submit the job
    start_dir = os.getcwd()
    os.chdir(irc_dir)
    irc_job = job_manager.SlurmJob()
    slurm_cmd = f"sbatch {slurm_run_file}"
    irc_job.submit(slurm_cmd)

    # only wait once all jobs have been submitted
    os.chdir(start_dir)
    irc_job.wait(check_interval=60)
