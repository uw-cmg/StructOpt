import os
import logging
import numpy as np
import shutil
import time
#from ovito.io import import_file
#from ovito.modifiers import CoordinationAnalysisModifier
#from asap3.analysis.rdf import RadialDistributionFunction
#from RadicalDistributionFunction_function import RadialDistributionFunction
from diffpy.structure import Structure
from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
from diffpy.structure.parsers import getParser

from ase.io import write
import gparameters
from structopt.io import write_xyz
from structopt.tools import root, single_core, parallel
from structopt.common.crossmodule.exceptions import FEMSIMError


class RDF(object):
    """Contains parameters and functions for running RDF through Python."""

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = self.read_inputs(parameters)
        # These parameteres do not need to exist between generations
        # They are used for before/after rdf processing
        self.base = None
        self.folder = None
        self.paramfilename = None

    @single_core
    def read_inputs(self, parameters):
        data = open(parameters.kwargs.rdf_data_filename).readlines()
        data.pop(0)  # Comment line
        data = [line.strip().split()[:3] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        r, gr= zip(*data)
        # Set r and gr data for chi2 comparison
        self.r = np.array(r)
        self.gr = np.array(gr)
        return parameters

    @single_core
#    def write_xyzfile(self, individual):
        # Write structure file to disk so that the fortran femsim can read it in
#        individual.set_cell([[self.parameters.kwargs.xsize, 0., 0.], [0., self.parameters.kwargs.ysize, 0.], [0., 0., self.parameters.kwargs.zsize]])
#        individual.wrap()
#        for index in range(0, 3):
#            lo = np.amin(individual.get_positions()[:, index])
#            hi = np.amax(individual.get_positions()[:, index])
#            assert lo >= 0
#            assert hi <= self.parameters.kwargs.xsize
#        comment = "{} {} {}".format(self.parameters.kwargs.xsize, self.parameters.kwargs.ysize, self.parameters.kwargs.zsize)
#        filename = os.path.join(gparameters.logging.path, 'modelfiles', 'individual{id}.xyz'.format(id=individual.id))
#        write_xyz(filename, individual, comment=comment)
       

    @single_core
    #rdf calculated by ovito
    #def get_gr_data(self, individual):
    #    filename = os.path.join(gparameters.logging.path, 'modelfiles', 'individual{id}.xyz'.format(id=individual.id))
    #    timeout = 10.  # seconds
    #    interval = 0.3  # seconds
    #    now = time.time()
    #    while True:
    #        if os.path.exists(filename):
    #            pipeline = import_file(filename)
    #            modifier = CoordinationAnalysisModifier(cutoff = self.r[-1], number_of_bins = len(self.r))
    #            pipeline.modifiers.append(modifier)
    #            data = pipeline.compute()
    #            gr = np.asarray(data.tables['coordination-rdf'].xy())[:, 1]
    #            break
    #        if time.time() - now > timeout:
    #            raise RDFError("rdf could not be read after trying for {} seconds.".format(timeout))
    #        time.sleep(interval)
    #    return gr
    def get_gr_data(self, individual):
        timeout = 300.  # seconds
        interval = 0.3  # seconds
        now = time.time()
        cfg = { 'qmax' : 30,
                'qdamp': 0.01,
                #'qstep': 0.01,
                'delta1': 0.002,
                'rmin' : 0.05,
                'rmax' : 10.05,
                'rstep': 0.05,
        }
        pc1 = DebyePDFCalculator(**cfg)
        tmp_cif = gparameters.logging.path + '/cifs'
        if not os.path.exists(tmp_cif):
            os.makedirs(tmp_cif)
        cif_file = tmp_cif+'/individual%s.cif'%individual.id
        individual.write(cif_file)
        idv_struct = Structure(filename=cif_file)
        while True:
            #gr = RadialDistributionFunction(atoms=individual, rMax=10.0, nBins=500).get_rdf(c="reduced")
            #gr = RadialDistributionFunction(atoms=individual, rMax=10.5, nBins=70).get_rdf(c="volume")   
            if self.parameters.kwargs.opt_rdflatt == True:
                
                latt_rdf = np.arange(self.parameters.kwargs.rdflatt_start, self.parameters.kwargs.rdflatt_end,  self.parameters.kwargs.rdflatt_step)
                fitness_rdf = []
                gr_lst = []
                for latt in latt_rdf:
                    idv_struct.lattice.a = idv_struct.lattice.b = idv_struct.lattice.c = latt
                    atom_num = len(idv_struct)
                    volume = idv_struct.lattice.volume
                    r1, gr_sim = pc1(idv_struct)
                    gr_sim = gr_sim-4*np.pi*r1*atom_num/volume
                    fitness_rdf_unit = individual.fitnesses.RDF.chi2(gr_sim)
                    fitness_rdf.append(fitness_rdf_unit)
                    gr_lst.append(gr_sim)
                min_idx = fitness_rdf.index(min(fitness_rdf))
                lattice_return =  latt_rdf[min_idx]
                fitness_return =  fitness_rdf[min_idx]
            else:
                atom_num = len(idv_struct)
                volume = idv_struct.lattice.volume
                r1, gr_sim = pc1(idv_struct)
                gr_sim = gr_sim-4*np.pi*r1*atom_num/volume
                fitness_rdf = individual.fitnesses.RDF.chi2(gr_sim)
                lattice_return =  self.parameters.kwargs.xsize
                fitness_return =  fitness_rdf
                
            break
        if time.time() - now > timeout:
            raise RDFError("rdf could not be read after trying for {} seconds.".format(timeout))
        time.sleep(interval)
        os.remove(cif_file)
        return lattice_return, fitness_return

    @single_core
    def chi2(self, gr):
        #return np.sum(((self.vk - vk) / self.vk_err)**2) / len(self.k)
        #return np.sum((self.gr - gr)**2) / len(self.r)
        sf = 1 #scaling factor 
        # fitness in r range
        rmax = 10
        nbin = 200
        r_ini = 1.65
        r_fin = 7.00
        start = int(r_ini/rmax*nbin) - 1
        end = int(r_fin/rmax*nbin)

        #exp_gr = self.gr
        #sim_gr = gr 

        #idx = np.round(np.linspace(0, len(gr[start:end])-1, 25)).astype(int)
        #exp_gr = self.gr[start:end][idx]
        #sim_gr = gr[start:end][idx]

        exp_gr = self.gr[start:end]
        sim_gr = gr[start:end]

        return np.sum((sim_gr/sf - exp_gr)**2) / len(sim_gr)
        

