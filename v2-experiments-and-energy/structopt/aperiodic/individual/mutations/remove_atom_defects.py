import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers

def remove_atom_defects(individual, surf_CN=11):
    """Moves atoms around on the surface based on coordination number
    Moves a surface atom with a low CN to an atom with a high CN

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    surf_CN : int
        The maximum coordination number to considered a surface atom

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if len(individual) == 0:
        return False

    # Analyze the individual
    CNs = CoordinationNumbers(individual)
    
    # Get indices, CNs, and positions of all surface sites
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs) if CN <= surf_CN]
    if len(surf_indices_CNs) == 0:
        return False
    surf_indices, surf_CNs = list(zip(*surf_indices_CNs))
    surf_indices = list(surf_indices)
    surf_CNs = list(surf_CNs)
    surf_positions = np.array([individual.positions[i] for i in surf_indices])
    surf_CN_counts = {CN: surf_CNs.count(CN) for CN in set(surf_CNs)}
    surf_probs = [2.0 ** -CN / surf_CN_counts[CN] for CN in surf_CNs]
    surf_probs /= sum(surf_probs)

    surf_index = np.random.choice(surf_indices, p=surf_probs)

    individual.pop(surf_index)

    return 
