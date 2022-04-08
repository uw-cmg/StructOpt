import logging

from structopt.tools import root, single_core, parallel
from structopt.tools.parallel import allgather
import gparameters


@parallel
def fitness(population, parameters):
    """Perform the LAMMPS fitness calculation on an entire population.

    Args:
        population (Population): the population to evaluate
    """
    rank = gparameters.mpi.rank
    if parameters.use_mpi4py:
        logger = logging.getLogger('by-rank')
    else:
        logger = logging.getLogger('output')

    to_fit = [individual for individual in population if not individual._fitted]

    if not to_fit:
        return [individual.LAMMPS for individual in population]

    if parameters.use_mpi4py:
        ncores = gparameters.mpi.ncores
    else:
        ncores = 1

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(to_fit):
        individuals_per_core[i % ncores].append(individual)

    for individual in individuals_per_core[rank]:
        print("Running LAMMPS fitness evaluation on individual {}".format(individual.id))
        energy = individual.fitnesses.LAMMPS.calculate_fitness(individual)
        convertor_kcal2ev = 0.043361254529175
        individual.LAMMPS = energy*convertor_kcal2ev
        logger.info('Individual {0} after LAMMPS evaluation has energy {1}'.format(individual.id, energy))

    fits = [individual.LAMMPS for individual in population]
    positions_per_core = {rank: [population.position(individual) for individual in individuals] for rank, individuals in individuals_per_core.items()}
    if parameters.use_mpi4py:
        fits = allgather(fits, positions_per_core)

    # Save the fitness value for the module to each individual after they have been allgathered
    for i, fit in enumerate(fits):
        population.get_by_position(i).LAMMPS = fit

    return [individual.LAMMPS for individual in population]

