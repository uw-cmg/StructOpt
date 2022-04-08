import numpy as np

from structopt.tools import root, single_core, parallel
from scipy.constants import physical_constants

@single_core
def roulette(fits, nkeep, T=None):
    """Select individuals with a probability proportional to their fitness.
    Fitnesses are renormalized from 0 - 1, which means minimum fitness
    individual is never included in in the new population.

    Parameters
    ----------
    fits : dict<int, float>
        Dictionary of <individual.id, fitness> pairs.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    T : float
        If T is not None, a boltzman-like transformation is applied
        to all fitness values with T.
    """

    ids, fits = zip(*fits.items())
    ids = list(ids)
    fits = list(fits)

    # Normalize fits from 0 (min fit) to 1 (max fit)
    fit_max = max(fits)
    fit_min = min(fits)
    best_id = ids[fits.index(fit_min)]

    if T is None:
        fits = np.array([-(fit - fit_max) for fit in fits])
        fits /= np.nan_to_num(max(fits))
    else:
        k = physical_constants['Boltzmann constant in eV/K'][0]
        fits = np.array([-(fit - fit_min) for fit in fits])
        fits = np.exp(fits/(k*T))

    # Generate probabilities and pick individuals
    p = np.nan_to_num(fits / np.sum(fits))

    # If less species have nonzero probability than nkeep, select
    # all nonzero probability and select random zero probability
    ids_nonzero_p = np.array([i for i, p_i in zip(ids, p) if p_i != 0])
    ids_zero_p = np.array([i for i, p_i in zip(ids, p) if p_i == 0])

    if len(ids_nonzero_p) < nkeep:
        n_zero_p_to_add = nkeep - len(ids_zero_p) - len(ids_nonzero_p)
        ids_zero_p_keep = np.random.choice(ids_zero_p, n_zero_p_to_add, replace=False)
        to_keep = np.append(ids_nonzero_p, ids_zero_p_keep)
    else:
        to_keep = np.random.choice(ids_nonzero_p, nkeep, replace=False, p=p)

    return to_keep

