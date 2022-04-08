import os
import logging
import numpy as np
import shutil
import time

import gparameters
from structopt.io import write_xyz
from structopt.tools import root, single_core, parallel
from structopt.common.crossmodule.exceptions import FEMSIMError


class FEMSIM(object):
    """Contains parameters and functions for running FEMSIM through Python."""

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = self.read_inputs(parameters)
        # Multiply the experimental data by the thickness scaling factor
        self.vk = np.multiply(self.parameters.kwargs.thickness_scaling_factor, self.vk)
        self.vk_err = np.multiply(self.parameters.kwargs.thickness_scaling_factor, self.vk_err)
        self.parameters.path = gparameters.logging.path

        # These parameteres do not need to exist between generations
        # They are used for before/after femsim processing
        self.base = None
        self.folder = None
        self.paramfilename = None

        assert self.parameters.kwargs.xsize == self.parameters.kwargs.ysize == self.parameters.kwargs.zsize


    @single_core
    def read_inputs(self, parameters):
        data = open(parameters.kwargs.vk_data_filename).readlines()
        data.pop(0)  # Comment line
        data = [line.strip().split()[:3] for line in data]
        data = [[float(line[0]), float(line[1]), float(line[2])] for line in data]
        k, vk, vk_err = zip(*data)
        # Set k and vk data for chi2 comparison
        self.k = np.array(k)
        self.vk = np.array(vk)
        self.vk_err = np.array(vk_err)
        return parameters


    @single_core
    def update_parameters(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self.parameters, key, value)


    @single_core
    def get_spawn_args(self, individual):
        """Returns a dictionary of arguments to be passed to MPI.COMM_SELF.Spawn which will be collected for all
        structures and concatenated into MPI.COMM_SELF.Spawn_multiple:
        https://github.com/mpi4py/mpi4py/blob/2acfc552c42846628304e54a3b87e2bf3a59af07/src/mpi4py/MPI/Comm.pyx#L1555
        """
        femsim_command = os.environ['FEMSIM_COMMAND']
        args = [self.base, self.paramfilename]
        info = {'wdir': self.folder}
        return {'command': femsim_command, 'args': args, 'info': info}


    @single_core
    def setup_individual_evaluation(self, individual):

        logger = logging.getLogger('by-rank')

        logger.info('Received individual HI = {0} for FEMSIM evaluation'.format(individual.id))

        # Make individual folder and copy files there
        # print('\n==========================')
        # print('FEMSIM generation {}, individual {}'.format(gparameters.generation, individual.id))
        # print('==========================\n')
        self.folder = os.path.abspath(os.path.join(self.parameters.path, 'FEMSIM/generation{gen}/individual{i}'.format(gen=gparameters.generation, i=individual.id)))
        os.makedirs(self.folder, exist_ok=True)
        if not os.path.isfile(os.path.join(self.folder, self.parameters.kwargs.vk_data_filename)):
            shutil.copy(self.parameters.kwargs.vk_data_filename, os.path.join(self.folder, self.parameters.kwargs.vk_data_filename))

        self.paramfilename = os.path.join(self.folder, "femsim.{}.in".format(individual.id))
        shutil.copy(self.parameters.kwargs.parameter_filename, self.paramfilename)
        self.write_paramfile(individual)

        self.base = 'indiv{i}'.format(i=individual.id)


    @single_core
    def write_paramfile(self, individual):
        # Write structure file to disk so that the fortran femsim can read it in
        individual.set_cell([[self.parameters.kwargs.xsize, 0., 0.], [0., self.parameters.kwargs.ysize, 0.], [0., 0., self.parameters.kwargs.zsize]])
        individual.wrap()
        for index in range(0, 3):
            lo = np.amin(individual.get_positions()[:, index])
            hi = np.amax(individual.get_positions()[:, index])
            #assert lo >= 0
            #assert hi <= self.parameters.kwargs.xsize
            if lo >= 0 or hi <= self.parameters.kwargs.xsize:
                print('coordination exceed the box!!! carefull with the FEMSIM simulation!!!')
        comment = "{} {} {}".format(self.parameters.kwargs.xsize, self.parameters.kwargs.ysize, self.parameters.kwargs.zsize)
        filename = os.path.join(gparameters.logging.path, 'modelfiles', 'individual{id}.xyz'.format(id=individual.id))
        write_xyz(filename, individual, comment=comment)

        with open(self.paramfilename, 'w') as f:
            f.write('# Parameter file for generation {gen}, individual {i}\n'.format(gen=gparameters.generation, i=individual.id))
            f.write('{}\n'.format(filename))
            f.write('{}\n'.format(self.parameters.kwargs.vk_data_filename))
            f.write('{}\n'.format(self.parameters.kwargs.Q))
            f.write('{} {} {}\n'.format(self.parameters.kwargs.nphi, self.parameters.kwargs.npsi, self.parameters.kwargs.ntheta))


    @single_core
    def get_vk_data(self):
        filename = os.path.join(self.folder, 'vk_initial_{base}.txt'.format(base=self.base))
        # filename = os.path.join('./', 'vk_initial_{base}.txt'.format(base=self.base))
        timeout = 10.  # seconds
        interval = 0.3  # seconds
        now = time.time()
        while True:
            if os.path.exists(filename):
                data = open(filename).readlines()
                if len(self.vk) == len(data) and data[-1][-1] == '\n':  # Then the file loaded with all the data
                    data = [line.strip().split()[:2] for line in data]
                    data = [[float(line[0]), float(line[1])] for line in data]
                    vk = np.array([vk for k, vk in data])
                    break
            if time.time() - now > timeout:
                raise FEMSIMError("V(k) could not be read after trying for {} seconds.".format(timeout))
            time.sleep(interval)
        return vk


    @single_core
    def chi2(self, vk):
        if self.parameters.kwargs.opt_beta == True:
            chi2_lst = []
            beta_lst = np.arange(self.parameters.kwargs.beta_start, self.parameters.kwargs.beta_end, self.parameters.kwargs.beta_step)
            for beta in beta_lst:
                vk_exp = self.vk * beta / self.parameters.kwargs.thickness_scaling_factor
                vk_err_exp = self.vk_err * beta / self.parameters.kwargs.thickness_scaling_factor
                chi2_lst.append(np.sum(((vk_exp - vk) / vk_err_exp)**2) / len(self.k))
            
            min_idx_femsim = chi2_lst.index(min(chi2_lst))
            return chi2_lst[min_idx_femsim], beta_lst[min_idx_femsim]
        else:
            return np.sum(((self.vk - vk) / self.vk_err)**2) / len(self.k), self.parameters.kwargs.thickness_scaling_factor
#        return np.sum((self.vk - vk)**2 / self.vk) / len(self.k)

