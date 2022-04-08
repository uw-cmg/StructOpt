from structopt.tools import root, single_core, parallel
import gparameters

import multiprocessing as mp
from multiprocessing import Process
import time

gb_proc_pool = None

def create_process_pool(njobs):
    global gb_proc_pool
    if gb_proc_pool is None:
        gb_proc_pool = mp.Pool(njobs)

def individual_runner(in_params):
    individual, gparams = in_params
    individual.relaxations.hard_sphere_cutoff.relax(individual)
    return individual

@parallel
def relax(population, parameters):
    """Relax the entire population using a hard-sphere cutoff method.

    Args:
        population (Population): the population to relax
    """

    njobs = parameters['njobs'] if 'njobs' in parameters else 1
    if njobs > 1:
        create_process_pool(njobs)

    to_relax = [individual for individual in population if not individual._relaxed]
    ncores = gparameters.mpi.ncores
    rank = gparameters.mpi.rank

    individuals_per_core = {rank: [] for rank in range(ncores)}
    for i, individual in enumerate(to_relax):
        individuals_per_core[i % ncores].append(individual)

    t0 = time.time()
    if njobs > 1:
        # When using mp.Pool.map, the input variables are copy by value 
        # in contrast with noraml copy by reference python function call.
        # So here we have to return the modified individual.
        # ref. https://stackoverflow.com/questions/57533533/using-multiprocessing-pool-map-to-change-the-value-of-variables-passed-to-it
        param_lst = [(indiv, gparameters) for indiv in to_relax]
        relaxed_lst = gb_proc_pool.map(individual_runner, param_lst)
        population.update(relaxed_lst)
        individuals_per_core[0] = relaxed_lst
    else:
        for individual in individuals_per_core[rank]:
            individual.relaxations.hard_sphere_cutoff.relax(individual)
    t1 = time.time()
    print('******* TOTAL HARD SPHERE OPTIMIZATION elapsed: {:.6f}s *******'.format(t1 - t0))
    
    if parameters.use_mpi4py:
        population.allgather(individuals_per_core)

