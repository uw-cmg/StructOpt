import scipy
import numpy as np
from copy import deepcopy


def tournament(population, fits, tournament_size=5, unique_pairs=True,
               unique_parents=True, keep_best=False):
    """Selects pairs in seperate "tournaments", where a subset of the
    population are randomly selected and the highest fitness allowed to pass.
    In addition to a population, their fits, and end population size, takes in
    a tournament size parameter.

    Parameters
    ----------
    population : Population
        The population of individuals needed to be trimmed
    fits : list
        List of fitnesses that correspond to the population.
    tournament_size : int
        The number of individuals in each tournament. If 1,
        tournament is the same as random selection. If
        len(population), corresponds to the "best" selection process
    unique_pairs : bool
        If True, all combinations of parents are unique, though parents
        can show up in different pairs. True increases the diversity of 
        the population.
    unique_parents : bool
        If True, all parents can only mate with on other individual.
        True increases the diversity of the population.
    """

    # Get ranks of each population value based on its fitness
    ranks = list(scipy.stats.rankdata(fits, method='ordinal'))
    ids_population = [individual.id for individual in population]
    id_to_rank = {id: rank for id, rank in zip(ids_population, ranks)}
    rank_to_id = {rank: id for id, rank in zip(ids_population, ranks)}
    n_pairs = int(len(fits) / 2)
    pairs_id = []

    for _ in range(n_pairs):
        # Perform tournament for father selection
        if len(ids_population) > tournament_size:
            tournament_father = np.random.choice(ids_population, tournament_size, replace=False)
        else:
            tournament_father = ids_population

        if keep_best:
            max_rank = min([id_to_rank[id] for id in ids_population])
        else:
            max_rank = min([id_to_rank[id] for id in tournament_father])
        id_father = rank_to_id[max_rank]

        # Choose the second parent based on renormalized probabilities
        # First remove the father from the population and make a temporary
        # population for choosing the mother. If unique_parents
        # is on, remove him from the global population as well.
        ids_population_temp = deepcopy(ids_population)
        ind_delete = ids_population.index(id_father)
        del ids_population_temp[ind_delete]

        if unique_parents:
            del ids_population[ind_delete]

        # Now remove mothers that would make repeat father/mother pairs
        if unique_pairs:
            for pair in deepcopy(pairs_id):
                if id_father in pair:
                    del pair[pair.index(id_father)]
                    id_mother = pair[0]
                    ind_delete = ids_population_temp.index(id_mother)
                    del ids_population_temp[ind_delete]

        # Using subset, perform tournament to find mother
        if len(ids_population_temp) > tournament_size:
            tournament_mother = np.random.choice(ids_population_temp, tournament_size, replace=False)
        else:
            tournament_mother = ids_population_temp
        max_rank = min([id_to_rank[id] for id in tournament_mother])
        id_mother = rank_to_id[max_rank]

        if unique_parents:
            ind_delete = ids_population.index(id_mother)
            del ids_population[ind_delete]

        pairs_id.append([id_father, id_mother])

    pairs = [[population[i], population[j]] for i, j in pairs_id]

    return pairs

