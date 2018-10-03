import logging

from structopt.tools import root, single_core, parallel
from structopt.tools.parallel import allgather
import gparameters

@parallel
def fitness(population, parameters):
    """Performs the STEM fitness calculation on an entire population.

    Parameters
    ----------
        population : structopt.Population
    """

    to_fit = [individual for individual in population if not individual._fitted]

    if parameters.use_mpi4py:
        logger = logging.getLogger('by-rank')
        ncores = gparameters.mpi.ncores
    else:
        logger = logging.getLogger('output')
        ncores = 1

    rank = gparameters.mpi.rank

    individuals_per_core = {r: [] for r in range(ncores)}
    for i, individual in enumerate(to_fit):
        individuals_per_core[i % ncores].append(individual)

    for individual in individuals_per_core[rank]:
        print("Evaluating fitness of individual {} on rank {} with STEM".format(individual.id, rank))
        chi2 = individual.fitnesses.STEM.calculate_fitness(individual)
        individual.STEM = chi2
        logger.info('Individual {0} after STEM evaluation has chi^2 {1}'.format(individual.id, chi2))

    positions_per_core = {rank: [population.position(individual) for individual in individuals] for rank, individuals in individuals_per_core.items()}

    fits = []
    for individual in population:
        if hasattr(individual, 'STEM'):
            fits.append(individual.STEM)
        else:
            fits.append(None)

    if parameters.use_mpi4py:
        fits = allgather(fits, positions_per_core)

    # Save the fitness value for the module to each individual after they have been allgathered
    for i, fit in enumerate(fits):
        population.get_by_position(i).STEM = fit

    return [individual.STEM for individual in population]    

