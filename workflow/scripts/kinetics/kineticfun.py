# Functions for running a kinetics job using this workflow
import pandas as pd
import re
import os
import glob
import subprocess
import autotst.reaction
import autotst.calculator.gaussian
import job_manager
from hotbit import Hotbit  # TODO - move this and other autoTST dependencies elsewhere
import ase.io.gaussian
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
        return f'{list_of_indices[0]}, {list_of_indices[1]}'

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

    reactant_token = reaction_smiles.split('_')[0]
    product_token = reaction_smiles.split('_')[1]

    reactant_tokens = reactant_token.split('+')
    product_tokens = product_token.split('+')

    # print(product_tokens)
    for reactant_str in reactant_tokens:
        if reactant_str == 'carbonmonoxide':
            reactant_str = '[C-]#[O+]'
        reactant = rmgpy.species.Species(smiles=reactant_str)
        reactants.append(reactant)
    for product_str in product_tokens:
        if product_str == 'carbonmonoxide':
            product_str = '[C-]#[O+]'
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
    1 for Error termination - due to all degrees of freedom being frozen
    2 for Error termination not falling under 1
    -1 for no termination
    """
    error_termination = False
    with open(log_file, 'rb') as f:
        f.seek(0, os.SEEK_END)
        normal_termination = False
        error_termination = False
        for i in range(0, 10):
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
                return 1
        if error_termination:
            return 2
        return -1


def run_TS_shell_calc(reaction_index, use_reverse=False):
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
            if status == 0 or status == 1:
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
    reaction.generate_conformers(ase_calculator=Hotbit())
    print('Done generating conformers in AutoTST...')
    print(f'{len(reaction.ts[direction])} conformers found')
    with open(logfile, 'a') as f:
        f.write('Done generating conformers in AutoTST...\n')
        f.write(f'{len(reaction.ts[direction])} conformers found' + '\n')

    # Do the shell calculations
    # write Gaussian input files
    slurm_array_idx = []
    for i in range(0, len(reaction.ts[direction])):
        if i not in incomplete_indices and len(shell_gaussian_logs) > 0:
            with open(logfile, 'a') as f:
                f.write(f'skipping completed shell {i}' + '\n')
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
    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        f'fname="{shell_label[:-4]}.com"' + '\n\n',

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


def run_TS_overall_calc(reaction_index, use_reverse=False):
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
    overall_gaussian_logs = glob.glob(os.path.join(overall_dir, overall_label[:-8] + '*.log'))
    incomplete_indices = []
    for overall_ts_log in overall_gaussian_logs:
        if os.path.exists(overall_ts_log):
            status = termination_status(overall_ts_log)
            if status == 0:
                print(f'Overall TS optimization already ran for reaction {reaction_index}')
                with open(logfile, 'a') as f:
                    f.write(f'Overall TS optimization already ran for reaction {reaction_index}' + '\n')
            else:
                matches = re.search(overall_label[:-8] + '([0-9]{4})', overall_ts_log)
                run_index = int(matches[1])
                incomplete_indices.append(run_index)
    if not incomplete_indices and len(overall_gaussian_logs) > 0:
        return True

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

    slurm_array_idx = []
    for i in range(0, len(reaction.ts[direction])):
        if i not in incomplete_indices and len(overall_gaussian_logs) > 0:
            print(f'Skipping completed index {i}')
            continue

        overall_label = overall_label[:-8] + f'{i:04}.log'
        shell_opt = os.path.join(shell_dir, overall_label)
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
        'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        'RUN_i=$(printf "%04.0f" $(($SLURM_ARRAY_TASK_ID)))\n',
        f'fname="{overall_label[:-4]}.com"' + '\n\n',

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
    #arkane_job.submit(slurm_cmd)

