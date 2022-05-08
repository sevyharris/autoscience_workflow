# Functions for running a kinetics job using this workflow
import pandas as pd
import os
import autotst.reaction
import autotst.calculator.gaussian
import job_manager
from hotbit import Hotbit
import ase.io.gaussian


try:
    DFT_DIR = os.environ['DFT_DIR']
except KeyError:
    DFT_DIR = '/work/westgroup/harris.se/autoscience/autoscience_workflow/results/dft'


def get_num_reactions():
    """Function to lookup number of reactions in the reaction_list.csv
    """
    reaction_csv = os.path.join(DFT_DIR, '..', '..', 'resources', 'reaction_list.csv')
    reaction_df = pd.read_csv(reaction_csv)
    return reaction_df.i.values[-1]


def reaction_index2smiles(reaction_index):
    """Function to return reactuib smiles given a reactuib index
    looks up the results in the reactuib_list.csv
    """
    reaction_csv = os.path.join(DFT_DIR, '..', '..', 'resources', 'reaction_list.csv')
    reaction_df = pd.read_csv(reaction_csv)
    reaction_smiles = reaction_df.SMILES.values[reaction_index]
    return reaction_smiles


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

    # Check for already finished shell logfile
    shell_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        shell_label = 'rev_ts_0000.log'
        direction = 'reverse'
    shell_opt_log = os.path.join(shell_dir, shell_label)

    if os.path.exists(shell_opt_log):
        status = termination_status(shell_opt_log)
        if status == 0 or status == 1:
            print('Shell optimization already ran')
            with open(logfile, 'a') as f:
                f.write('Shell optimization already ran\n')
            return True

    reaction_smiles = reaction_index2smiles(reaction_index)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    reaction.ts[direction][0].get_molecules()
    reaction.generate_conformers(ase_calculator=Hotbit())

    # Do the shell calculation
    # write Gaussian input files
    if len(reaction.ts[direction]) > 1:
        raise ValueError(f'multiple {direction} ts')

    ts = reaction.ts[direction][0]
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
    slurm_run_file = os.path.join(shell_dir, 'run.sh')
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
        # TODO make this an array if multiple forward ts's
    }
    slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
    slurm_file_writer.settings = slurm_settings
    slurm_file_writer.content = [
        'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
        'mkdir -p $GAUSS_SCRDIR\n',
        'module load gaussian/g16\n',
        'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

        f'g16 {shell_label[:-4]}.com' + '\n',
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
    os.chdir(start_dir)
    shell_job.wait(check_interval=600)

    # TODO restart optimization if it fails. for now, move on


def run_TS_overall_calc(reaction_index, use_reverse=False):
    """Start a TS optimization from the geometry of the shell calculation
    """
    reaction_base_dir = os.path.join(DFT_DIR, 'kinetics', f'reaction_{reaction_index:04}')
    os.makedirs(reaction_base_dir, exist_ok=True)

    shell_dir = os.path.join(reaction_base_dir, 'shell')
    overall_dir = os.path.join(reaction_base_dir, 'overall')
    os.makedirs(overall_dir, exist_ok=True)

    overall_label = 'fwd_ts_0000.log'
    direction = 'forward'
    if use_reverse:
        overall_label = 'rev_ts_0000.log'
        direction = 'reverse'
    overall_ts_log = os.path.join(overall_dir, overall_label)

    # check that the overall job hasn't finished
    if os.path.exists(overall_ts_log):
        status = termination_status(overall_ts_log)
        if status == 0:
            print('Overall TS optimization already ran')
            with open(logfile, 'a') as f:
                f.write('Overall TS optimization already ran\n')
            return True

    reaction_smiles = reaction_index2smiles(reaction_index)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    reaction.ts[direction][0].get_molecules()
    reaction.generate_conformers(ase_calculator=Hotbit())

    if len(reaction.ts[direction]) > 1:
        raise ValueError(f'multiple {direction} ts')

    shell_opt = os.path.join(shell_dir, overall_label)
    try:
        with open(shell_opt, 'r') as f:
            atoms = ase.io.gaussian.read_gaussian_out(f)
            reaction.ts[direction][0]._ase_molecule = atoms
    except IndexError:
        # handle case where all degrees of freedom were frozen in the shell calculation
        if len(reaction.ts[direction][0]._ase_molecule) > 3:
            raise ValueError('Shell optimization failed to converge. Rerun it!')

    for i, ts in enumerate(reaction.ts[direction]):
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
        slurm_run_file = os.path.join(overall_dir, 'run.sh')
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
            # TODO make this an array if multiple forward ts's
        }

        slurm_file_writer = job_manager.SlurmJobFile(full_path=slurm_run_file)
        slurm_file_writer.settings = slurm_settings
        slurm_file_writer.content = [
            'export GAUSS_SCRDIR=/scratch/harris.se/guassian_scratch\n',
            'mkdir -p $GAUSS_SCRDIR\n',
            'module load gaussian/g16\n',
            'source /shared/centos7/gaussian/g16/bsd/g16.profile\n\n',

            f'g16 {overall_label[:-4]}.com' + '\n',
        ]
        slurm_file_writer.write_file()

        # submit the job
        start_dir = os.getcwd()
        os.chdir(overall_dir)
        overall_job = job_manager.SlurmJob()
        slurm_cmd = f"sbatch {slurm_run_file}"
        overall_job.submit(slurm_cmd)
        os.chdir(start_dir)
        overall_job.wait(check_interval=600)
