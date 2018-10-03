import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers

def remove_atom_random(individual, surf_CN=11):
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
    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN]
    surf_index = np.random.choice(surf_indices)

    individual.pop(surf_index)

    return 
