from structopt.tools import root, single_core, parallel
import gparameters


@parallel
def relax(population, parameters):
    """Relax the entire population using LAMMPS.

    Args:
        population (Population): the population to relax
    """

    to_relax = [individual for individual in population if not individual._relaxed]
    ncores = gparameters.mpi.ncores
    rank = gparameters.mpi.rank

    individuals_per_core = {rank: [] for rank in range(ncores)}
    for i, individual in enumerate(to_relax):
        individuals_per_core[i % ncores].append(individual)

    for individual in individuals_per_core[rank]:
        individual.relaxations.LAMMPS.relax(individual)

    if parameters.use_mpi4py:
        population.allgather(individuals_per_core)

