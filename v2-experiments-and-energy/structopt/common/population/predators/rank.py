import numpy as np
import scipy


def rank(fits, nkeep, p_min=None):
    """Selection function that chooses pairs of structures
    based on linear ranking.

    See "Grefenstette and Baker 1989 Whitley 1989".

    Parameters
    ----------
    fits : dict<int, float>
        Dictionary of <individual.id, fitness> pairs.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    p_min : float
        The probability of choosing the lowest ranked individual.
        Given population of size N, this should be below 1/nindiv.
        The probability of selecting rank N (worst) to rank 1 (best)
        increases from p_min to (2/N - p_min) in even, (1/N - p_min)
        increments. Defaults to (1/N)^2.
    """

    if p_min is None:
        p_min = 1.0 / len(fits) ** 2

    # Get ranks of each individual value based on its fitness
    ids, fits = zip(*fits.items())
    ids = list(ids)
    fits = list(fits)
    ranks = scipy.stats.rankdata(fits, method='ordinal')

    # Get probabilities based on linear ranking
    N = len(fits)
    eta_min = p_min * N
    eta_max = 2 - eta_min
    p_max = eta_max / N
    p = p_min + (p_max - p_min)*(N - ranks)/(N - 1)  # `ranks` is a numpy array, so p is too

    # Randomly choose `nkeep` values from the list `ids` given probability `p` for each value in `ids` (no duplicates)
    return  np.random.choice(ids, nkeep, replace=False, p=p)

