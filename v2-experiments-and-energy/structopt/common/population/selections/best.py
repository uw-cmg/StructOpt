import scipy.stats


def best(population, fits):
    """Deterministic selection function that chooses adjacently
    ranked individuals as pairs.

    Parameters
    ----------
    population : Population
        An population of individuals
    fits : list
        Fitnesses that corresponds to population
    """

    # Get ranks of each population value based on its fitness
    ranks = list(scipy.stats.rankdata(fits, method='ordinal'))
    ids = [individual.id for individual in population]
    pairs = []

    for i in range(len(population) // 2):
        rank_father = 2*i + 1
        rank_mother = 2*i + 2
        id_father = ids[ranks.index(rank_father)]
        id_mother = ids[ranks.index(rank_mother)]
        father = population[id_father]
        mother = population[id_mother]
        pairs.append([father, mother])

    return pairs

