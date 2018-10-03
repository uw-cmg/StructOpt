import scipy
import numpy as np


def tournament(fits, nkeep, tournament_size=5):
    """Selects individuals in seperate "tournaments", where a subset of the
    population are randomly selected and the highest fitness allowed to pass.
    In addition to a population, their fits, and end population size, takes in
    a tournament size parameter.

    Parameters
    ----------
    fits : dict<int, float>
        Dictionary of <individual.id, fitness> pairs.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    tournament_size : int
        The number of individuals in each tournament. If 1,
        tournament is the same as random selection. If
        len(population), corresponds to the "best" selection process
    """

    # Implementation NOT taken from:  https://en.wikipedia.org/wiki/Tournament_selection -- actually, this isn't the process we use.
    #   I also don't know how they define "p" in that pseudocode, but the variable `p_min` in rank.py might give a clue.
    # Implementation taken from:  Genetic Algorithms, Tournament Selection, and the Effects of Noise (http://www.complex-systems.com/pdf/09-3-2.pdf)
    # Another reference:  A Comparison of Selection Schemes Used in Evolutionary Algorithms (http://www.tik.ee.ethz.ch/file/6c0e384dceb283cd4301339a895b72b8/TIK-Report11.pdf)

    fits = fits.copy()

    to_keep = []
    for _ in range(nkeep):
        # choose k (the tournament size) individuals from the population at random
        if len(fits) > tournament_size:
            tournament_ids = np.random.choice(list(fits.keys()),
                                              size=tournament_size,
                                              replace=False)
        else:
            tournament_ids = list(fits.keys())

        # From a list of fitnesses from the tournament pool, get the id with the lowest fitness
        best_id = tournament_ids[ np.argmin([fits[id] for id in tournament_ids]) ]
        # Remove that id from the dictionary so we don't get it again in a later iteration
        fits.pop(best_id)
        to_keep.append(best_id)

    return to_keep

