import os
try:
    from ase import Atom, Atoms
    from ase.optimize import BFGS
    from ase.units import GPa
    from ase.calculators.neighborlist import NeighborList
except ImportError:
    pass
from StructOpt.inp_out.write_xyz import write_xyz
from StructOpt.fingerprinting import get_fingerprint
from StructOpt.tools.lammps import LAMMPS
from StructOpt.tools.setup_energy_calculator import setup_energy_calculator
from StructOpt.tools.compose_structure import compose_structure
from StructOpt.tools.decompose_structure import decompose_structure
import numpy
import math
import json
try:
    from mpi4py import MPI
except ImportError:
    pass
import logging
import pdb
import shutil
import time
import scipy
import random
#from StructOpt.tools.parallel_mpi4py import parallel_mpi4py
#from StructOpt.tools.serial_for import serial_for

class LAMMPS_eval(object):
    def __init__(self):
        
        self.args = self.read_inputs()
        for k,v in self.args.items():
            setattr(self,k,v)

            
    def read_inputs(self):
        args = dict()
        #args['parallel'] = True
        args = json.load(open('lammps_inp.json'))
        return args

    def evaluate_fitness(self, Optimizer, individ, relax=False):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        
        if rank==0:
            ntimes=int(math.ceil(1.*len(individ)/comm.Get_size()))
            
            nadd=int(ntimes*comm.Get_size()-len(individ))
            maplist=[[] for n in range(ntimes)]
            strt=0
            for i in range(len(maplist)):
                maplist[i]=[indi for indi in individ[strt:comm.Get_size()+strt]]
                strt+=comm.Get_size()
            for i in range(nadd):
                maplist[len(maplist)-1].append(None)
        else:
            ntimes=None
        ntimes = comm.bcast(ntimes,root=0)
        outs=[]

        for i in range(ntimes):
            if rank==0:
                one=maplist[i]
            else:
                one=None
            ind = comm.scatter(one,root=0)

            if ind == None:
                rank = MPI.COMM_WORLD.Get_rank()
                stro = 'Evaluated none individual on {0}\n'.format(rank)
                out = (None, stro)
            else:
                out = self.evaluate_indiv(Optimizer, ind, rank, relax)

            outt = comm.gather(out,root=0)
            if rank == 0:
                outs.extend(outt)
        return outs
    
    def evaluate_indiv(self, Optimizer, individ, rank, relax):

        logger = logging.getLogger(Optimizer.loggername)
        if relax:
            logger.info('Received individual HI = {0} for LAMMPS structure relaxation'.format(
                individ.history_index))
        else:
            logger.info('Received individual HI = {0} for LAMMPS energy evaluation'.format(
            individ.history_index))
        
        STR='----Individual ' + str(individ.history_index)+ ' Optimization----\n'
        indiv=individ[0]
        if 'EE' in Optimizer.debug:
            debug = True
        else:
            debug = False
        if debug: 
            write_xyz(Optimizer.debugfile,indiv,'Received by evaluate_fitness')
            Optimizer.debugfile.flush()
            logger.debug('Writing received individual to debug file')
        
        
        totalsol = compose_structure(Optimizer,individ)

        # Set calculator to use to get forces/energies
        if Optimizer.parallel:
            calc = self.setup_lammps(Optimizer, self.args, relax)
            if Optimizer.fixed_region:
                if debug:
                    logger.info('Setting up fixed region calculator')
                pms=copy.deepcopy(calc.parameters)
                try:
                    pms['mass'][len(pms['mass'])-1] += '\ngroup RO id >= {0}\nfix freeze RO setforce 0.0 0.0 0.0\n'.format(nat)
                except KeyError:
                    pms['pair_coeff'][0] += '\ngroup RO id >= {0}\nfix freeze RO setforce 0.0 0.0 0.0\n'.format(nat)
                calc = LAMMPS(parameters=pms, files=calc.files, keep_tmp_files=calc.keep_tmp_files, tmp_dir=calc.tmp_dir)
                lmin = copy.copy(self.args['minimize'])
                if debug:
                    logger.info('Setting up no local minimization calculator')
                self.args['minimize'] = None
                self.args['static_calc'] = self.setup_lammps(Optimizer, self.args, relax)
                self.args['minimize'] = lmin
        else:
            calc=Optimizer.calc
        totalsol.set_calculator(calc)
        totalsol.set_pbc(True)
        
        # Perform Energy Minimization
        if not Optimizer.parallel:
            if debug: 
                write_xyz(Optimizer.debugfile,totalsol,'Individual sent to Energy Minimizer')
                logger.debug('Writing structure sent to energy minimizer')
        try:
            cwd = os.getcwd()
            if debug:
                logger.info('Running local energy calculator')
            if Optimizer.fixed_region:
                totalsol, pea, energy, pressure, volume, STR = self.run_energy_eval(totalsol, Optimizer.fixed_region, STR, self.args['static_calc'])
            else:
                totalsol, pea, energy, pressure, volume, STR = self.run_energy_eval(totalsol, False, STR)
                logger.info('M:finish run_energy_eval, energy = {0} @ rank ={1}'.format(energy,rank))
        except Exception, e:
            logger.critical('Error in energy evaluation: {0}'.format(e), exc_info=True)
            path = os.path.join(cwd,'TroubledLammps')
            if not os.path.exists(path):
                os.mkdir(path)
            #Copy files over
            shutil.copyfile(calc.trajfile,os.path.join(path,os.path.basename(calc.trajfile)))
            shutil.copyfile(calc.infile,os.path.join(path,os.path.basename(calc.infile)))
            shutil.copyfile(calc.logfile,os.path.join(path,os.path.basename(calc.logfile)))
            shutil.copyfile(calc.datafile,os.path.join(path,os.path.basename(calc.datafile)))
            raise RuntimeError('{0}:{1}'.format(Exception,e))
        if not Optimizer.parallel:
            if debug:
                write_xyz(Optimizer.debugfile,totalsol,'Individual after Energy Minimization')
                Optimizer.debugfile.flush()
                logger.debug('Writing structure received from energy minimizer')
       
        
        individ, bul = decompose_structure(Optimizer,totalsol,individ)
        
        # Add concentration energy dependence
        if Optimizer.forcing=='energy_bias':
            if debug:
                logger.info('Applying energy bias for atoms with different number of atoms of type than in atomlist')
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                n[i]=len([inds for inds in totalsol if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=abs(n[i]-Optimizer.atomlist[i][1])
            factor=sum(n)**3
            energy=(energy+factor)/totalsol.get_number_of_atoms()
            STR+='Energy with Bias = {0}\n'.format(energy)
        elif Optimizer.forcing=='chem_pot':
            if debug:
                logger.info('Applying chemical potential bias for atoms with different number of atoms of type than in atomlist')
            n=[0]*len(Optimizer.atomlist)
            for i in range(len(Optimizer.atomlist)):
                n[i]=len([inds for inds in totalsol if inds.symbol==Optimizer.atomlist[i][0]])
                n[i]=n[i]*Optimizer.atomlist[i][3]
            factor=sum(n)
            energy=(energy+factor)/totalsol.get_number_of_atoms()
            STR+='Energy with Chemical Potential = {0}\n'.format(energy)

        individ.buli = bul
        individ.energy = energy
        individ.pressure = pressure
        individ.volume = volume

        ##Add pealist to include atom index based on sorted PE. 
        #logger.info('before sort{0}'.format(individ.energy))
        #self.sort_pealist(Optimizer,individ,pea)
        #energy = individ.energy
        #logger.info('after sort {0}'.format(individ.energy)) 
        if Optimizer.fingerprinting:
            if debug:
                logger.info('Identifying fingerprint of new structure')
            individ.fingerprint=get_fingerprint(Optimizer,individ,Optimizer.fpbin,Optimizer.fpcutoff)
        
        calc.clean()
        if relax:
            signal = 'Relaxed structure of individual {0} on {1}\n'.format(individ.index,rank)
        else:
            signal = 'Evaluated fitness of individual {0} on {1}\n'.format(individ.index,rank)
        signal += STR

        #self.fitness = energy
        if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
            individ.bulki = bul
       

        if relax:
            return individ, signal
        else:
            return individ.energy, signal

    def setup_lammps(self, Optimizer, args, relax):
        return setup_energy_calculator(Optimizer,'LAMMPS',relax)
        

    def update_parameters(self, **kwargs):
        # TODO
        #self.args['pair_coeff'] = 'eam'
        for key, value in kwargs.items():
            self.args[key] = value


    def sort_pealist(self, Optimizer,individ,pea):
        logger = logging.getLogger(Optimizer.loggername)

        #Add pealist to include atom index based on sorted PE. 
        syms = [sym for sym,c,m,u in Optimizer.atomlist]
        numatom = [c for sym,c,m,u in Optimizer.atomlist]
        peatom = [u for sym,c,m,u in Optimizer.atomlist]

        hpealist = []
        lpealist = []
        for j in range(len(syms)) :
            sym = syms[j]
            pelist = []
            peindexlist = []
            for i in range(len(pea)) :
              if individ[0][i].symbol == sym:
                if pea[i][0] < peatom[j] - 5.0 or pea[i][0] > peatom[j] + 5.0 :
                   individ.energy = 10000
                   message = 'Warning: Found oddly large energy from Lammps in structure HI={0}'.format(individ.history_index)
                   logger.warn(message)
                pelist.append(pea[i][0])
                peindexlist.append(i)
            pearray = numpy.asarray(pelist)
            pearray_sorted = pearray.argsort()
            for i in range(1,11):
               hpealist.append(peindexlist[pearray_sorted[-i]])
            for i in range(0,10):
               lpealist.append(peindexlist[pearray_sorted[i]])

            individ.hpealist = hpealist
            individ.lpealist = lpealist
        return
    

    def run_energy_eval(self, totalsol, fx_region=False, STR='', static_calc=None):
        totcop = totalsol.copy()
        OUT = totalsol.calc.calculate(totalsol)
        totalsol = OUT['atoms']
        pea = OUT['pea']
        totalsol.set_pbc(True)
        if fx_region:
            STR+='Energy of fixed region calc = {0}\n'.format(OUT['thermo'][-1]['pe'])
            totalsol.set_calculator(static_calc)
            OUT=totalsol.calc.calculate(totalsol)
            totalsol=OUT['atoms']
            totalsol.set_pbc(True)
            STR+='Energy of static calc = {0}\n'.format(OUT['thermo'][-1]['pe'])
        en=OUT['thermo'][-1]['pe']
        stress=numpy.array([OUT['thermo'][-1][i] for i in ('pxx','pyy','pzz','pyz','pxz','pxy')])*(-1e-4*GPa)
        volume = totalsol.get_volume()
        pressure = 0 ## should be modified if enthalpy_fit
        energy=en
        STR+='Energy per atom = {0}\n'.format(energy/len(totalsol))
        return totalsol, pea, energy, pressure, volume, STR
