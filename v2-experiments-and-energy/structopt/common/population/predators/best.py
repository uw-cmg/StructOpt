
def best(fits, nkeep):
    """Sorts individuals by fitness and keeps the top nkeep fitnesses.
    
    Parameters
    ----------
    fits : dict<int, float>
        Dictionary of <individual.id, fitness> pairs.
    nkeep : int
        The number of individuals to keep. In a GA run, corresponds
        to the sum of each generators number_of_individuals
    """
    sorted_ids, sorted_fits = zip(*sorted(fits.items(), key=lambda pair: pair[1]))
    return sorted_ids[:nkeep]

