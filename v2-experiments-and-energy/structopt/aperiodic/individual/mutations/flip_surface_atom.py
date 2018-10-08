import random
import numpy as np

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii

from ase.io import write

def flip_surface_atom(individual, surf_CN=11, cutoff=0.5):
    """Randomly "flips" an atom from one side of the particle to the other side
    of the particle or to a defect within a column.

    Parameters
    ----------
    individual : structopt.Individual object
        The individual object to be modified in place
    surf_CN : int
        The maximum coordination number for determining an atom as a surface atom
    cutoff : float
        The factor of the average bond length for searching in the xy directions for
        atoms in the same column

    Output
    ------
    out : None
        Modifies individual in-place
    """

    if len(individual) == 0:
        return False

    avg_bond_length = get_avg_radii(individual) * 2
    cutoff = avg_bond_length * cutoff
    
    # Analyze the individual
    CNs = CoordinationNumbers(individual)
    
    # Get all surface atoms
    positions = individual.get_positions()
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs)
                        if CN <= surf_CN and CN > 2]
    surf_indices, surf_CNs = list(zip(*surf_indices_CNs))
    surf_indices = list(surf_indices)
    surf_CNs = list(surf_CNs)

    # Pair each surface atom with a surface atom on the other side

    # First calculate the xy and z distances of each surface atom
    # with the other surface atoms. The xy coordinates are needed
    # to see if they're in the same column. The z coordinates are
    # needed to see if it is the top and or bottom atom in the column
    surf_positions = np.array([positions[i] for i in surf_indices])    
    surf_xys = np.array([surf_positions[:,:2]])
    surf_xy_vecs = surf_xys - np.transpose(surf_xys, [1, 0, 2])
    surf_xy_dists = np.linalg.norm(surf_xy_vecs, axis=2)

    surf_zs = np.array([surf_positions[:,-1]])
    surf_z_vecs = surf_zs - surf_zs.T
    surf_z_dists = np.absolute(surf_z_vecs)

    # For columns with multiple or zero surface atoms, we need to find
    # which one is on the other side of the particle
    flips = [np.where(dists < cutoff) for dists in surf_xy_dists]
    for i, column in enumerate(flips):
        column_z_vecs = surf_z_vecs[i][column]
        if (len(column_z_vecs) == 1
            or (not all(column_z_vecs >= 0) and not all(column_z_vecs <= 0))):
            flips[i] = None
        else:
            flips[i] = surf_indices[column[0][np.argmax(surf_z_dists[i][column])]]

    # Update the surface atom arrays
    for i, flip in reversed(list(enumerate(flips))):
        if flip is None:
            surf_positions = np.delete(surf_positions, i, axis=0)
            del surf_indices[i]
            del surf_CNs[i]
            del flips[i]

    flip_indices = [[surf_indices[i], j] for i, j in enumerate(flips)]
    flip_CN_diffs = [CNs[j] - CNs[i] for i, j in flip_indices]

    # Simple rank based selection based on coordination differences. Simply put
    # a move is more likely if a surface atom with a low coordination number
    # is moved on top of an atom with a high coordination number.
    min_CN_diff = np.min(flip_CN_diffs) - 1
    flip_CN_diffs -= min_CN_diff
    unique_CN_diffs = np.array(list(set(flip_CN_diffs)))
    flip_CN_ps =  2.0 ** unique_CN_diffs
    flip_CN_ps /= np.sum(flip_CN_ps)
    flip_CN = np.random.choice(unique_CN_diffs, p=flip_CN_ps)
    
    flip_indices_is = np.where(flip_CN_diffs == flip_CN)[0]
    flip_indices_i = random.choice(flip_indices_is)
    move_index, surf_index = flip_indices[flip_indices_i]

    # Add the atom at move_index to ontop/bottom of the surf_index
    move_atom = individual[move_index]
    surf_atom = individual[surf_index]
    if move_atom.z > surf_atom.z:
        move_atom.z = surf_atom.z - avg_bond_length
    else:
        move_atom.z = surf_atom.z + avg_bond_length

    return
