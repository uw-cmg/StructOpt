import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers

def enrich_surface(individual, surf_CN=11, species=None):
    """Mutation that selectively enriches the surface with a species.

    Parameters
    ----------
    individual : Individual
        An individual
    surf_CN : int
        The maximum coordination number of an atom to be considered surface
    species : str
        The surface to enrich with. If None, takes the lowest concentration
        species
    """

    syms = individual.get_chemical_symbols()
    if species is None:
        unique_syms = list(set(syms))
        counts = [syms.count(sym) for sym in unique_syms]
        species = unique_syms[np.argmin(counts)]

    # Get a random surface site that is not the species and a bulk site that is
    CNs = CoordinationNumbers(individual)

    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN and syms[i] != species]
    bulk_indices = [i for i, CN in enumerate(CNs) if CN > surf_CN and syms[i] == species]

    if len(surf_indices) == 0 or len(bulk_indices) == 0:
        return False

    surf_index = random.choice(surf_indices)
    surf_symbol = syms[surf_index]    
    bulk_index = random.choice(bulk_indices)

    individual[surf_index].symbol = species
    individual[bulk_index].symbol = surf_symbol    
    
    return
