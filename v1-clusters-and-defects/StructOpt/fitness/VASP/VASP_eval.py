from MAST.controllers.mastinput import MASTInput
from StructOpt.tools.setup_energy_calculator import setup_energy_calculator
import os
import glob
import shutil
from subprocess import Popen, PIPE, check_output
import json
from ase.io.vasp import write_vasp, read_vasp
try:
    from mpi4py import MPI
except ImportError:
    pass
from StructOpt.tools.compose_structure import compose_structure
from StructOpt.tools.decompose_structure import decompose_structure
import logging
import time

class VASP_eval(object):
    def __init__(self):
        self.args = self.read_inputs()
        for k,v in self.args.items():
            setattr(self,k,v)

    def read_inputs(self):
        args = json.load(open('vasp_inp.json'))
        return args

    def setup_mast_inp(self, Optimizer, individ, args, relax):
        
        tmp = os.path.join(os.environ["STRUCTOPT_INSTALL_PATH"],"StructOpt/fitness/VASP/mast_template.inp")
        with open(tmp,'r') as fp:
            mastfile = fp.read()
        
        mf = mastfile.split('\n')
        for l in range(len(mf)):
            if 'system_name' in mf[l]:
                mf[l] += ' ('
                for ind in range(len(individ)):
                    if ind > 0:
                        mf[l] += ', '
                    mf[l] += 'Indiv%s'%ind
                mf[l] += ')'

            if (relax and ('begin relax'==mf[l])) or ((not relax) and ('begin static'==mf[l])):
                for k,v in args.items():
                    mf[l] += '\n%s  %s'%(k,v) 
            if '$recipe' in mf[l]:
                mf[l+1] += ' (%s)'%('relax' if relax else 'static')
            if 'posfile' in mf[l]:
                mf[l] += ' ('
                for ind in range(len(individ)):
                    if ind > 0:
                        mf[l] += ', '
                    mf[l] += 'POSCAR_Indiv%s'%ind
                mf[l] += ')'
 
        with open('my_mast.inp','w') as f:
            f.write('\n'.join(mf))
        return

    def setup_vasp(self, Optimizer, args):
        return setup_energy_calculator(Optimizer,'VASP')

    def evaluate_fitness(self, Optimizer, individ, relax=False):
        rank = MPI.COMM_WORLD.Get_rank()

        if rank == 0:
            cwd = os.getcwd()
            try:
                os.mkdir('{filename}-rank0/VASPFiles'.format(filename=Optimizer.filename))
            except OSError:
                pass
            os.chdir('{filename}-rank0/VASPFiles'.format(filename=Optimizer.filename))
            out = self.evaluate_indiv(Optimizer, individ, relax)
            out = zip(*out)
            os.chdir(cwd)
        else:
            out = None
        out = MPI.COMM_WORLD.bcast(out, root=0)
        return out
    
    def evaluate_indiv(self, Optimizer, individ, relax):
        logger = logging.getLogger(Optimizer.loggername)
        logger.info('Setting up MAST input for individual evaluation with VASP')

        self.setup_mast_inp(Optimizer, individ, self.args, relax)
        num = 0
        for ind in individ:
            totalsol = compose_structure(Optimizer,ind)
            write_vasp('POSCAR_Indiv%s'%num, totalsol, direct=True, vasp5=True)
            num += 1
        #try:
        #    os.mkdir('scratch')
        #    os.mkdir('archive')
        #except OSError:
        #    pass

        cwd = os.getcwd()
        #scratch = os.path.join(cwd,'scratch')
        #archive = os.path.join(cwd,'archive')
        
        #os.environ['MAST_SCRATCH'] = scratch
        #os.environ['MAST_ARCHIVE'] = archive
        scratch = os.getenv('MAST_SCRATCH')
        archive = os.getenv('MAST_ARCHIVE')

        #proc = Popen(['mast','-i','my_mast.inp'])
        #proc.wait()

        mast = MASTInput(inputfile="my_mast.inp")
        ind_folders = mast.check_independent_loops()
        

        for files in glob.glob("loop_*.inp"):
            os.remove(files)
        
        logger.info('Completed setting up MAST jobs for individual evaluation of generation #%s'%Optimizer.generation)
        logger.info("Type 'mast' to manually submit VASP jobs or wait for crontab for submission.")

        while True:

            #proc = Popen(['mast'])
            #proc.wait()
            stats = 'VASP jobs stats:'
            finished = 'Finished jobs: '
            complete = 0
            for d in ind_folders:
                if d in os.listdir(archive):
                    finished += '%s, '%d
                    complete += 1
            running = 'Running jobs: '
            for d in ind_folders:
                if d in os.listdir(scratch):
                    running += '%s, '%d
            logger.info(stats)
            logger.info(finished)
            logger.info(running)
            if complete == len(individ):
               logger.info('Finished VASP jobs for all individuals of generation #%s'%Optimizer.generation) 
               break
            time.sleep(60)
        
        stro = []
        energies = []
        for i in range(len(ind_folders)): 
            ingred = os.path.join(archive,ind_folders[i],'Individual')
            energy_output = check_output(['tail','-n1',os.path.join(ingred, 'OSZICAR')])
            individ[i].energy = float(energy_output.split()[4])
            energies.append(individ[i].energy)
            if relax:
                stro.append('Relaxed structure of individual %s\n'%individ[i].history_index)
            else:
                stro.append('Evaluated energy of individual %s to be %s\n'%(individ[i].history_index,individ[i].energy))
            totalsol = read_vasp(os.path.join(ingred, 'CONTCAR'))
            individ[i], bul = decompose_structure(Optimizer,totalsol,individ[i])
        
            shutil.rmtree(os.path.join(archive,ind_folders[i])) # either remove or zip ??

        

        if relax:
            return individ, stro
        return energies, stro
