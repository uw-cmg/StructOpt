import random
import numpy as np
from ase import Atoms

from structopt.common.individual import Individual
from structopt.tools import root, single_core, parallel
from structopt.tools import random_three_vector
from structopt.common.crossmodule import repair_cluster
from structopt.common.crossmodule.similarity import get_offset

import time

@single_core
def rotate_fixed(individual1, individual2, align_cutoff=0.1, repair_composition=True):
    """Similar to rotate except the children aren't borne out of cut out
    versions.

    Parameters
    ----------
    individual1 : Individual
        The first parent
    individual2 : Individual 
        The second parent
    center_at_atom : bool
        This centers the cut at an atom. This is particularly useful 
        when one desires a crystalline solution. If both parents
        are crystalline, the children will likely not have grain boundaries.
    repair_composition : bool 
        If True, conserves composition. For crossovers that create children
        with more (less) atoms, atoms are taken from (added to) the surface
        of the particle. For incorrect atomic ratios, atomic symbols are
        randomly interchanged throughout the particle

    Returns
    -------
    Individual: The first child
    Individual: The second child
    """

    # Translate individuals so atoms are most aligned
    

    # Pick a random rotation angle and vector
    rot_vec = random_three_vector()
    rot_angle = random.uniform(0, 180) * np.pi / 180.0

    # Rotate both individuals
    ind1c.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))
    ind2c.rotate(rot_vec, a=rot_angle, center=(0, 0, 0))

    # Create the children
    child1 = Atoms(cell=ind1c.get_cell())
    child2 = Atoms(cell=ind2c.get_cell())
    child1.extend(ind1c[ind1c.positions[:,2] >= 0])
    child1.extend(ind2c[ind2c.positions[:,2] < 0])
    child2.extend(ind2c[ind2c.positions[:,2] >= 0])
    child2.extend(ind1c[ind1c.positions[:,2] < 0])

    # Repair the children
    if repair_composition:
        syms = ind1c.get_chemical_symbols()
        atomlist = [[sym, syms.count(sym)] for sym in set(syms)]
        repair_cluster(child1, atomlist)
        repair_cluster(child2, atomlist)

    # Reorient the children
    child1.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child1.center()
    child2.rotate(rot_vec, a=-rot_angle, center=(0, 0, 0))
    child2.center()

    # Convert the children to an Individual with the parent
    # module parameters
    full_child1 = ind1c.copy(include_atoms=False)
    full_child1.extend(child1)
    full_child2 = ind2c.copy(include_atoms=False)
    full_child2.extend(child2)
    
    return full_child1, full_child2
