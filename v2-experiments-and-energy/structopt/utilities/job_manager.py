"""This contains a StructOpt calculator for submitting jobs to queue, tracking their progress, and reading their output."""

import os
import re
import json
import shutil
from copy import deepcopy
from operator import itemgetter
import subprocess
from subprocess import PIPE
from importlib import import_module

import numpy as np
import ase

import structopt.utilities
from ..common.individual import Individual
from .exceptions import UnknownState, Running, Queued, Submitted
from .data_explorer.core import DataExplorer

class JobManager(object):
    '''This class is responsible for the submission and tracking of jobs'''

    def __init__(self, calcdir=None, optimizer='genetic.py', parameters=None,
                 submit_parameters={'system':'PBS'}):
        if parameters is None:
            parameters = {}

        # Initialize inputs
        if calcdir is None:
            self.calcdir = os.getcwd()
        else:
            self.calcdir = os.path.expandvars(calcdir)

        self.log_dir = None

        if not os.path.isfile(os.path.abspath(optimizer)):
            optimizer = os.path.expandvars('$STRUCTOPT_HOME/structopt/optimizers/{}'.format(optimizer))
        else:
            optimizer = os.path.abspath(os.path.expandvars(optimizer))

        self.optimizer = optimizer
        self.parameters = deepcopy(parameters)
        self.submit_parameters = deepcopy(submit_parameters)
        if 'job_name' not in self.submit_parameters:
            self.submit_parameters['job_name'] = self.calcdir

        self.path = os.path.abspath(self.calcdir)
        self.cwd = os.getcwd()
        self.system_name = os.path.basename(self.calcdir)

        # If totally clean calculation
        if not os.path.isdir(self.path):
            self.status = 'clean'
            os.makedirs(self.path)

        # If directory exists
        elif not os.path.isfile(os.path.join(self.path, 'structopt.in.json')):
            self.status = 'clean'

        # If input exists but was never run. This also applies to cancelled
        # jobs that were never submitted
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and not self.job_in_queue(os.path.join(self.path, 'jobid'))
              and not self.read_runs()):
            self.read_input()
            self.status = 'initialized'

        # If the job is running
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and self.job_in_queue(os.path.join(self.path, 'jobid'))):
            self.read_input()
            if self.status == 'running':
                if not self.read_runs():
                    self.status = 'queued'
                else:
                    self.read_generations()

        # If the job is done, check the output
        elif (os.path.isfile(os.path.join(self.path, 'structopt.in.json'))
              and not self.job_in_queue(os.path.join(self.path, 'jobid'))
              and self.read_runs()):
            self.read_input()
            self.read_generations()
            if self.generations is not None:
                self.status = 'done'
            else:
                self.status = 'error'

        self.parameters.update(parameters)

    def restart(self): # TODO
        """Loads up the last generation of a previous run and modifies 
        the self.parameters to load up those structures on the next run"""

        XYZs_dir = os.path.join(self.log_dir, 'XYZs/generation{}'.format(self.generations[-1]))
        fnames = [os.path.join(XYZs_dir, f) for f in os.listdir(XYZs_dir) if f.endswith('.xyz')]
        new_generator = {'generators': {'read_extxyz': {'number_of_individuals': len(fnames),
                                                        'kwargs': fnames}}}
        self.parameters.update(new_generator)
        self.status = 'initialized'

        return

    def get_jobid(self, jobid=None):
        if jobid is None:
            jobid = os.path.join(self.path, 'jobid')

        if not os.path.exists(jobid):
            return False

        with open(jobid) as f:
            jobid = f.readline().split()[-1]

        return jobid

    def job_in_queue(self, jobid='jobid'):
        '''return True or False if the directory has a job in the queue'''
        jobid = self.get_jobid(jobid)
        if jobid is False:
            return False

        # Behavior will depend on whether we are in the slurm or pbs environment
        if self.submit_parameters['system'] == 'PBS':
            try:
                jobids_in_queue = subprocess.check_output('qselect')
            except FileNotFoundError:
                self.status = 'done'
                return False
            jobids_in_queue = [job.decode('utf-8') for job in jobids_in_queue.split()]
        else:
            raise NotImplemented(self.submit_parameters['system'], 'not implemented yet')

        if jobid in jobids_in_queue:
            # get details on specific jobid
            output, error = subprocess.Popen(['qstat', format(jobid)], stdout=PIPE).communicate()
            if error is None:
                fields = output.decode('utf-8').split('\n')[-2].split()
                job_status = fields[4]
                if job_status == 'C':
                    return False
                elif job_status == 'R':
                    self.status = 'running'
                    return True
                else:
                    self.status = 'queued'
                    return True
            return False
        else:
            return False


    ####################################################################
    ### Calculation methods. Includes write, run, and submit scripts ###
    ####################################################################

    def optimize(self, run_method='submit', rerun=False, restart=False):
        """Runs the optimizer"""

        run_method = getattr(self, run_method)

        if restart and self.status == 'done':
            self.restart()

        if self.status in ['clean', 'initialized'] or rerun:
            run_method()
        elif self.status == 'running':
            raise Running
        elif self.status == 'queued':
            raise Queued

        return

    def run(self): # TODO
        """Runs the job as is in the current directory."""

        self.write_input()
        self.write_submit()

        return

    def submit(self):
        """Submits the job to the queue. Do this in the calculation
        directory"""

        from .rc import QUEUE_OPTIONS as queue_options

        self.write_input()
        shutil.copy(self.optimizer, os.path.join(self.path, os.path.basename(self.optimizer)))
        self.write_submit()

        submit_cmd = queue_options[self.submit_parameters['system']]['submit']

        os.chdir(self.path)

        p = subprocess.Popen([submit_cmd, 'submit.sh'], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()

        with open('jobid', 'wb') as f:
            f.write(out)

        os.chdir(self.cwd)

        raise Submitted(out)

    def write_input(self):
        """Writes the parameters input to a json file"""

        input_file = os.path.join(self.path, 'structopt.in.json')
        with open(input_file, 'w') as f:
            json.dump(self.parameters, f, indent=4, sort_keys=True)

        return

    def write_submit(self):
        """This function writes the submit script."""

        from .rc import QUEUE_OPTIONS as queue_options
        from .rc import RUN_OPTIONS as run_options
        from .rc import CUSTOM_LINES as custom_lines

        # Gather variables from rc file
        submit = self.submit_parameters
        queue_system = submit['system']
        options = queue_options[queue_system]

        prefix = options['prefix']
        job_name = options['job_name'].format(submit['job_name'])
        queue = options['queue'].format(submit['queue'])
        nodes_cores = options['nodes_cores'].format(submit['nodes'], submit['cores'])
        total_cores = submit['nodes'] * submit['cores']
        walltime = options['walltime'].format(submit['walltime'])
        misc = options['misc']

        mpirun = run_options['mpirun']
        python = run_options['python']
        optimizer = os.path.basename(self.optimizer)
        input_file = 'structopt.in.json'

        # Write the submit script
        script = '#!/bin/bash\n\n'

        script += '{prefix} {job_name}\n'.format(**locals())
        script += '{prefix} {queue}\n'.format(**locals())
        script += '{prefix} {nodes_cores}\n'.format(**locals())
        script += '{prefix} {walltime}\n'.format(**locals())
        script += '{prefix} {misc}\n\n'.format(**locals())

        script += '{custom_lines}\n\n'.format(**locals())

        script += 'cd {}\n\n'.format(self.path)

        if total_cores > 1:
            script += '{mpirun} -n {total_cores} {python} {optimizer} {input_file}'.format(**locals())
        else:
            script += '{python} {optimizer} {input_file}'.format(**locals())

        submit_file = os.path.join(self.path, 'submit.sh')
        with open(submit_file, 'w') as f:
            f.write(script)

        return

    def read_runs(self):
        """Stores the output directory in the order in which they were run. 
        Sets StructOpt to read results form most recent by default"""

        dirs = [d for d in os.listdir(self.path)
                if os.path.isdir(os.path.join(self.path, d))]
        pattern = r'logs(.*)'
        log_times = [int(re.match(pattern, d, re.I|re.M).group(1)) for d in dirs
                     if re.match(pattern, d, re.I|re.M)]
        log_dirs = [os.path.join(self.path, 'logs{}'.format(t)) for t in log_times]

        if len(log_dirs) > 0:
            log_dirs_times = zip(log_dirs, log_times)
            log_dirs_times = sorted(log_dirs_times, key=itemgetter(1))
            log_dirs, log_times = zip(*log_dirs_times)
            log_dirs = list(log_dirs)
            log_times = list(log_times)
            for i, log_dir in reversed(list(enumerate(log_dirs))):
                if 'fitnesses.log' not in os.listdir(log_dirs[i]):
                    log_dirs.pop(i)
            self.log_dirs = log_dirs
            self.set_run(-1)
            return True
        else:
            self.log_dirs = log_dirs
            return False

    def get_number_of_runs(self):
        return len(self.log_dirs)

    def set_run(self, run_number):
        """Sets the get and read functions on a certain run number.

        Parameters
        ----------
        run_number : int
            The run number we wish to extra data out of. Normal indexing rules
            apply, so run_number = -1 is the last (most recent) run, while 0
            is the first (earliest) run.
        """

        if self.log_dirs is None:
            self.read_runs()

        new_log_dir = self.log_dirs[run_number]

        self.log_dir = new_log_dir
        self.read_generations()

    def check_run(self):
        """Check the stdout to see if a run is complete. NOTE, only works
        for jobs submitted to the queue and jobs where the stdout is
        saved to stdout.txt"""

        if 'stdout.txt' in os.listdir(self.path):
            out_file = os.path.join(self.path, 'stdout.txt')
        else:
            files = os.listdir(self.path)
            pattern = '.*.o(.*)'
            files = [f for f in files if re.match(pattern, f, re.I|re.M)]
            files = [f for f in files if re.match(pattern, f, re.I|re.M).group(1).isnumeric()]
            if len(files) == 0:
                return 'error'
            jobids = [int(re.match(pattern, f, re.I|re.M).group(1)) for f in files]
            files_jobids = zip(files, jobids)
            files_jobids = sorted(files_jobids, key=lambda i: i[1])
            out_file = files_jobids[-1][0]

        with open(os.path.join(self.path, out_file)) as f:
            last = None
            for last in (line for line in f if line.rstrip('\n')):
                pass

        # Check if the job ran error free. Happens with "Finished!" gets
        # printed and if the job was killed due to walltime
        if ('Finished!' in last):
            return 'done'
        elif ('walltime' in last):
            return 'timeout'
        else:
            return 'error'

    def read_generations(self):
        """Determines which generations exists and initializes
        a list to store the individuals"""

        # Get list of generations available to be read for output
        if not os.path.isdir(os.path.join(self.log_dir, 'XYZs')):
            return 
        dirs = [d for d in os.listdir(os.path.join(self.log_dir, 'XYZs'))]
        pattern = r'generation(.*)'
        generations = [int(re.match(pattern, d, re.I|re.M).group(1)) for d in dirs]
        self.generations = sorted(generations)

        # Initialize dictionary for storing XYZ coordinates
        self.populations = [None for generation in range(max(generations) + 1)]

        return

    def read_input(self): # TODO
        """Read the input and store them self.parameter"""

        with open(os.path.join(self.path, 'structopt.in.json')) as f:
            parameters = json.load(f)

        self.parameters.update(parameters)

        return

    def get_data_explorer(self, run_number=-1):
        """Returns a DataExplorer instance for analysis of a specific run"""

        log_dir = self.log_dirs[run_number]
        return DataExplorer(log_dir)
