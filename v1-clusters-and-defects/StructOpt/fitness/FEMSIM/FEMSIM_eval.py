import subprocess
import shlex
import sys
import json
import time
import numpy as np
import os
from StructOpt.inp_out.write_xyz import write_xyz
import logging
import math
import shutil
from mpi4py import MPI

class FEMSIM_eval(object):
    def __init__(self):
        self.args = self.read_inputs()

        #self.step_number = 0

        self.vk = np.multiply(self.args['thickness_scaling_factor'], self.vk)  # Multiply the experimental data by the thickness scaling factor


    def read_inputs(self):
        args = json.load(open('femsim_inp.json'))

        data = open(args['vk_data_filename']).readlines()
        data.pop(0)  # Comment line
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        k, vk = zip(*data)
        # Set k and vk data for chi2 comparison
        self.k = np.array(k)
        self.vk = np.array(vk)
        return args


    def update_parameters(self, **kwargs):
        #self.step_number += 1

        #self.args['aberation_coef'] += 1

        for key, value in kwargs.items():
            self.args[key] = value

    def evaluate_fitness(self, Optimizer, individ):
        rank = MPI.COMM_WORLD.Get_rank()
        out = []
        if rank==0:
            femsimfiles = '{filename}-rank0/FEMSIMFiles'.format(filename=Optimizer.filename)
            try:
                os.mkdir(femsimfiles)
            except OSError:
                pass
            
            for i in range(len(individ)):
                indiv_folder = '{filename}-rank0/FEMSIMFiles/Individual{i}'.format(filename=Optimizer.filename,i=i)
                try:
                    os.mkdir(indiv_folder)
                except OSError:
                    pass
                if not os.path.isfile(os.path.join(indiv_folder,self.args['vk_data_filename'])):
                    shutil.copy(self.args['vk_data_filename'], os.path.join(indiv_folder,self.args['vk_data_filename']))
                cwd = os.getcwd()
                os.chdir(indiv_folder)
                out.append(self.evaluate_indiv(Optimizer, individ[i], i))
                os.chdir(cwd)

        out = MPI.COMM_WORLD.bcast(out, root=0)
        return out

    def evaluate_indiv(self, Optimizer, individ, i):

        logger = logging.getLogger(Optimizer.loggername)

        logger.info('Received individual HI = {0} for FEMSIM evaluation'.format(
            individ.history_index))
        
        
        paramfilename = self.args['parameter_filename']
        paramfilename = paramfilename.split('.')
        paramfilename[-2] = '{head}_{i}'.format(head=paramfilename[-2], i=individ.history_index) 
        paramfilename = '.'.join(paramfilename)

        self.write_paramfile(paramfilename, Optimizer, individ, i)

        base = 'indiv{i}'.format(i=individ.history_index)
        self.run_femsim(base=base, paramfilename=paramfilename)
        vk = self.get_vk_data(base)

        chisq = self.chi2(vk)
        logger.info('M:finish chi2 evaluation, chi2 = {0}'.format(chisq))
        stro = 'Evaluated individual {0}\n'.format(individ.history_index)
        
        return chisq, stro


    def write_paramfile(self, paramfilename, Optimizer, individ, i):
        # Write structure file to disk so that the fortran femsim can read it in
        #ase.io.write('structure_{i}.xyz'.format(i=individ.history_index), individ[0])
        data = "{} {} {}".format(self.args['xsize'], self.args['ysize'], self.args['zsize'])
        write_xyz('structure_{i}.xyz'.format(i=individ.history_index), individ[0], data)
        
        with open(paramfilename, 'w') as f:
            f.write('# Parameter file for generation {gen}, individual {i}\n'.format(gen=Optimizer.generation, i=individ.history_index))
            f.write('{}\n'.format('structure_{i}.xyz'.format(i=individ.history_index)))
            f.write('{}\n'.format(self.args['vk_data_filename']))
            f.write('{}\n'.format(self.args['Q']))
            f.write('{} {} {}\n'.format(self.args['nphi'], self.args['npsi'], self.args['ntheta']))
            f.write('{}\n'.format(self.args['thickness_scaling_factor']))



    def run_femsim(self, base, paramfilename):
        self.run_subproc('{femsim_command} {base} {paramfilename}'.format(femsim_command=os.getenv('FEMSIM_COMMAND'),base=base, paramfilename=paramfilename))


    def run_subproc(self, args):
        """ args should be the string that you would normally run from bash """
        #print("Running (via python): {0}".format(args))
        sargs = shlex.split(args)
        p = subprocess.Popen(sargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = []
        for nextline in iter(p.stdout.readline, ""):
            sys.stdout.write(nextline)
            output.append(nextline)
            sys.stdout.flush()
        poutput = p.stdout.read()
        perr = p.stderr.read()
        preturncode = p.wait()
        if(preturncode != 0):
            print("{0} exit status: {1}".format(args, preturncode))
            print("{0} failed: {1}".format(args, perr))
        return ''.join(output)

    def get_vk_data(self, base):
        # Sleep until we can get the file
        # There may be in issue if the file is only partially written to when the open command gets run...
        # Let's hope that doesn't happen. If it does, I will convert this to a `if f in os.listdir()` command, followed by another short sleep.
        while True:
            try:
                data = open('vk_initial_{base}.txt'.format(base=base)).readlines()
                break
            except IOError:
                time.sleep(5.0)
        #data.pop(0)  # I don't see the comment line in the vk output data
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk


    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)

