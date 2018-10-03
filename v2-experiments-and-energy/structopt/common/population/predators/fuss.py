import random
import numpy as np


def fuss(fits, nkeep, nbest=0, fusslimit=10):
    """Fixed uniform selection scheme. Aimed at maintaining diversity
    in the population. In the case where low fit is the highest
    fitness, selects a fitness between min(fits) and min(fits) + fusslimit,
    if the difference between the min(fit) and max(fit) is larger than fusslimit.

    Parameters
    ----------
    fits : dict<int, float>
        Dictionary of <individual.id, fitness> pairs.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    nbest : int
        The top n individuals to always keep (default 0)
    fusslimit : float
        Individuals that have fitness fusslimit
        worse than the max fitness will not be considered
    """

    # Find min and max fitness
    minf = min(fits.values())
    maxf = max(fits.values())
    if abs(maxf-minf) > fusslimit:
            maxf = minf + fusslimit

    # Select random point on fitness line
    pt = random.uniform(minf, maxf)

    to_keep = []
    # Always keep the top nbest individuals
    if nbest > 0:
        fits = sorted(fits.items(), key=lambda pair: pair[1])
        sorted_ids, sorted_fits = zip(*fits)
        to_keep.extend(sorted_ids[:nbest])

    # Calculate the distance of each individual's fitness from that point
    distances = {id: np.absolute(fit - pt) for id, fit in fits.items()}

    # Select individuals with lowest distance (ie closest to the selected point)
    distances = sorted(distances.items(), key=lambda pair: pair[1])
    sorted_ids, sorted_distances = zip(*distances)
    to_keep.extend(sorted_ids[:nkeep - nbest])

    return to_keep

