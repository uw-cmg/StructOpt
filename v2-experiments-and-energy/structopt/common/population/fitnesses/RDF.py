import logging
import math
from collections import defaultdict
import subprocess
import os
import shutil
import numpy as np

from structopt.tools.parallel import root, single_core, parallel, parse_MPMD_cores_per_structure
import gparameters

import multiprocessing as mp
from multiprocessing import Process
import time

# njobs set the maximum number of parallel relaxation jobs (aka femsim jobs).
# The number of actual parallel jobs will be min(#. of to_fit indiviuals, njobs)

gb_proc_pool = None

def create_process_pool(njobs):
    global gb_proc_pool
    if gb_proc_pool is None:
        gb_proc_pool = mp.Pool(njobs)

def rdf_runner(indiv):
    t0 = time.time()
    indiv.rdflatt, indiv.RDF = indiv.fitnesses.RDF.get_gr_data(indiv)
    t1 = time.time()
    print('******* Individual {} rdf elapsed: {:.6f}s *******'.format(indiv.id, t1 - t0))
    return indiv

@root
def fitness(population, parameters):
    """Perform the FEMSIM fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """

    njobs = parameters['njobs'] if 'njobs' in parameters else 1
    if njobs > 1:
        create_process_pool(njobs)

    from mpi4py import MPI    

    to_fit = [individual for individual in population if not individual._fitted]
    if parameters.skip_bad_lammps and all(hasattr(individual, "LAMMPS") for individual in population):
        to_fit = [individual for individual in to_fit if individual.LAMMPS != np.inf]

    t0 = time.time()
    while to_fit:
        if njobs > 1:
            fitted_indiv_lst = gb_proc_pool.map(rdf_runner, to_fit)
            population.update(fitted_indiv_lst)
            break
        else:
            for indiv in to_fit:
                indiv.rdflatt, indiv.RDF = indiv.fitnesses.RDF.get_gr_data(indiv)
            break
    t1 = time.time()
    print('******* TOTAL RDF elapsed: {:.6f}s *******'.format(t1 - t0))

    return [individual.RDF for individual in population]

