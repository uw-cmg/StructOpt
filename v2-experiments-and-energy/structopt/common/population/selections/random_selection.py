from itertools import combinations
import random


def random_selection(population, fits):
    """Randomly selects parents
    
    Parameters
    ----------
    population : Population
        An population of individuals
    fits : list
        Fitnesses that corresponds to population
    """

    pairs = [pair for pair in combinations(population, 2)]
    random.shuffle(pairs)
    pairs = pairs[:len(population) // 2]
    return pairs

