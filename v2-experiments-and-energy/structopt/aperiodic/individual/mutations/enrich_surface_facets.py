import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers

def enrich_surface_facets(individual, surf_CN=11, species=None):
    """Mutation that selectively enriches facets with a species. Facets
    are defined as atoms atoms with higher coordination numbers.

    Parameters
    ----------
    surf_CN : int
        The maximum coordination number of an atom to be considered surface
    species : str
        The surface to enrich with. If None, takes the lowest concentration
    """

    syms = individual.get_chemical_symbols()
    if species is None:
        unique_syms = list(set(syms))
        counts = [syms.count(sym) for sym in unique_syms]
        species = unique_syms[np.argmin(counts)]

    # Get a random surface site that is not the species and a bulk site that is
    CNs = CoordinationNumbers(individual)

    defect_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN and syms[i] == species]
    defect_CNs = [CNs[i] for i in defect_indices]
    facet_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN and syms[i] != species]
    facet_CNs = [CNs[i] for i in facet_indices]

    if len(defect_indices) == 0 or len(facet_indices) == 0:
        return False

    unique_defect_CNs = np.unique(defect_CNs)
    defect_probs = 2.0 ** (surf_CN + 1 - unique_defect_CNs)
    defect_probs /= np.sum(defect_probs)
    defect_CN = np.random.choice(unique_defect_CNs, p=defect_probs)
    defect_index = defect_indices[random.choice(np.where(np.array(defect_CNs) == defect_CN)[0])]

    unique_facet_CNs = np.unique(facet_CNs)
    facet_probs = 2.0 ** (unique_facet_CNs)
    facet_probs /= np.sum(facet_probs)
    facet_CN = np.random.choice(unique_facet_CNs, p=facet_probs)
    facet_index = facet_indices[random.choice(np.where(np.array(facet_CNs) == facet_CN)[0])]
    facet_symbol = syms[facet_index]

    individual[defect_index].symbol = facet_symbol
    individual[facet_index].symbol = species
    
    return
