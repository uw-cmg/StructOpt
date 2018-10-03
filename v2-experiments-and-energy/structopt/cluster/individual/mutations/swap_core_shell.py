import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers

def swap_core_shell(individual, surf_CN=11):
    """Swaps atoms on the surface with an atom in the core. Only does it
    for different element types.
    
    Parameters
    ----------
    individual : Individual
        An individual
    surf_CN : int
        The maximum coordination number of an atom to be considered surface
    """

    if not len(individual):
        return None

    CNs = CoordinationNumbers(individual)

    surf_indices = [i for i, CN in enumerate(CNs) if CN < surf_CN and CN > 3]
    bulk_indices = [i for i in range(len(individual)) if i not in surf_indices]

    # Construct surface and bulk dictionaries of elements and their indices
    syms = individual.get_chemical_symbols()
    surf_elements = list(set([syms[i] for i in surf_indices]))
    bulk_elements = list(set([syms[i] for i in bulk_indices]))

    surf_dict = {element: [] for element in surf_elements}
    bulk_dict = {element: [] for element in bulk_elements}

    for i in surf_indices:
        surf_dict[syms[i]].append(i)

    for i in bulk_indices:
        bulk_dict[syms[i]].append(i)

    # Get a list of bulk indices that CAN be swapped for each
    # unique surface element
    swap_dict = {}
    for surf_element in surf_dict:
        swap_dict[surf_element] = []
        for bulk_element in bulk_dict:
            if bulk_element == surf_element:
                continue
            swap_dict[surf_element] += bulk_dict[bulk_element]

    # Get a list of all swaps, taken based on their probability of happening        
    swap_list = []

    # First pick an element to swap based on the concentration in
    # the available surface swap atoms and atoms it can swap with
    surf_prob = np.asarray([len(surf_dict[element]) * len(swap_dict[element])
                            for element in surf_elements], dtype=float)
    if np.sum(surf_prob) == 0:
        return
    surf_prob /= np.sum(surf_prob)
    surf_element = np.random.choice(surf_elements, p=surf_prob)
    surf_index = np.random.choice(surf_dict[surf_element])

    # Pick a random bulk atom to swap to
    bulk_index = np.random.choice(swap_dict[surf_element])
    bulk_element = individual[bulk_index].symbol

    # Swap the elements
    individual[surf_index].symbol = bulk_element
    individual[bulk_index].symbol = surf_element

    return None
