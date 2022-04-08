import logging
import math
from collections import defaultdict
import subprocess
import os
import shutil

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

def femsim_runner(in_params):
    indiv_jobname, indiv_file = in_params
    t0 = time.time()
    femsim_command = os.environ['FEMSIM_COMMAND']
    subprocess.run([femsim_command, indiv_jobname, indiv_file])

    t1 = time.time()
    print('******* Individual {} femsim elapsed: {:.6f}s *******'.format(indiv_jobname, t1 - t0))
    return

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

    for individual in to_fit:
        individual.fitnesses.FEMSIM.setup_individual_evaluation(individual)

    t0 = time.time()
    while to_fit:
        if njobs > 1: 
            # use mp.Pool
            in_params = [(indiv.fitnesses.FEMSIM.base, indiv.fitnesses.FEMSIM.paramfilename) for indiv in to_fit]
            gb_proc_pool.map(femsim_runner, in_params)
            for indiv in to_fit:
                folder = indiv.fitnesses.FEMSIM.folder
                jobname = indiv.fitnesses.FEMSIM.base
                shutil.move(os.path.join('./', 'vk_initial_{base}.txt'.format(base=jobname)), \
                    os.path.join(folder, 'vk_initial_{base}.txt'.format(base=jobname)))
                vk = indiv.fitnesses.FEMSIM.get_vk_data()
                indiv.FEMSIM, indiv.beta = indiv.fitnesses.FEMSIM.chi2(vk)
            break
        ncores = gparameters.mpi.ncores
        cores_per_individual = ncores // len(to_fit)
        # Round cores_per_individual down to nearest power of 2
        if cores_per_individual == 0:
            # We have more individuals than cores, so each fitness scheme needs to be run multiple times
            cores_per_individual = 1

        pow(2.0, math.floor(math.log2(cores_per_individual)))
        minmax = parse_MPMD_cores_per_structure(parameters.MPMD)
        assert minmax['max'] >= minmax['min'] > 0
        if cores_per_individual < minmax['min']:
            cores_per_individual = minmax['min']
        elif cores_per_individual > minmax['max']:
            cores_per_individual = minmax['max']

        logger = logging.getLogger('output')

        # Setup each individual and get the inputs for each individual that need to be passed into the spawn
        multiple_spawn_args = defaultdict(list)
        for individual in to_fit:
            spawn_args = individual.fitnesses.FEMSIM.get_spawn_args(individual)
            for arg_name, arg in spawn_args.items():
                multiple_spawn_args[arg_name].append(arg)

        # Make sure each key in `multiple_spawn_args` has the same number of elements and set `count` to that value
        count = -1
        for key, value in multiple_spawn_args.items():
            if count != -1:
                assert len(value) == count
            count = len(value)

        # Create MPI.Info objects from the kwargs dicts in multiple_spawn_args['info'] for each rank
        infos = [MPI.Info.Create() for _ in multiple_spawn_args['info']]
        for i, info in enumerate(infos):
            for key, value in multiple_spawn_args['info'][i].items():
                info.Set(key, value)

        # Run the multiple spawn
        individuals_per_iteration = max(1, ncores // cores_per_individual)  # Error correction if the user allowed more cores for MPMD than they are running on
        individuals_per_iteration = min(individuals_per_iteration, len(to_fit))
        num_iterations = math.ceil(len(to_fit) / individuals_per_iteration)
        for i in range(num_iterations):
            j = i * individuals_per_iteration
            individuals_this_iteration = individuals_per_iteration
            if i == len(to_fit) // individuals_per_iteration:  # The last iteration may not be exactly individuals_per_iteration
                individuals_this_iteration = len(to_fit) % individuals_per_iteration
                cores_per_individual = ncores // individuals_this_iteration # All the cores available should be used
            print("Spawning {} femsim processes, each with {} cores".format(individuals_this_iteration, cores_per_individual))
            # intercomm = MPI.COMM_SELF.Spawn_multiple(command=multiple_spawn_args['command'][j:j+individuals_this_iteration],
            #                                          args=multiple_spawn_args['args'][j:j+individuals_this_iteration],
            #                                          maxprocs=[cores_per_individual]*individuals_this_iteration,
            #                                          info=infos[j:j+individuals_this_iteration]
            #                                          )
            # # Disconnect the child processes
            # intercomm.Disconnect()
            ncores = 1
            subprocess.run(['mpirun', '-np', str(ncores)] \
                + multiple_spawn_args['command'][j:j+individuals_this_iteration] \
                + multiple_spawn_args['args'][j:j+individuals_this_iteration][0])
            

            # Collect the results for each chisq and return them
            for i, individual in enumerate(to_fit[j:j+individuals_this_iteration]):
                shutil.move(os.path.join('./', 'vk_initial_{base}.txt'.format(base=individual.fitnesses.FEMSIM.base)), \
                    os.path.join(individual.fitnesses.FEMSIM.folder, 'vk_initial_{base}.txt'.format(base=individual.fitnesses.FEMSIM.base)))
                vk = individual.fitnesses.FEMSIM.get_vk_data()
                individual.FEMSIM = individual.fitnesses.FEMSIM.chi2(vk)
                logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, individual.FEMSIM))
        break
    t1 = time.time()
    print('******* TOTAL FEMSIM elapsed: {:.6f}s *******'.format(t1 - t0))
    # print([individual.FEMSIM for individual in population])

    return [individual.FEMSIM for individual in population]

