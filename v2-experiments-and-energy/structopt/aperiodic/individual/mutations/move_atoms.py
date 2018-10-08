import random
import numpy as np

from structopt.common.crossmodule import NeighborList

def move_atoms(individual, max_natoms=0.20):
    """Randomly moves atoms within a cluster.
    
    Parameters
    ----------
    individual : Individual 
        An individual
    max_natoms : float or int
        if float, the maximum number of atoms that will be moved is 
        max_natoms*len(individual). if int, the maximum number of atoms 
        that will be moved is max_natoms default: 0.20
    """

    if not len(individual):
        return False

    if isinstance(max_natoms, float):
        assert max_natoms <= 1
        max_natoms = int(len(individual)*max_natoms)

    NNs = NeighborList(individual)
    CNs = [len(NN) for NN in NNs]
    
    # Get the surface atoms. These provide bounds for the moves
    # First get unit vectors and mags of all surface atoms
    positions = individual.get_positions()
    com = np.sum(positions.T, axis=1) / len(individual)
    surf_indices = [i for i, CN in enumerate(CNs) if CN < 11]
    surf_positions = np.array([positions[i] for i in surf_indices])
    surf_magnitudes = np.linalg.norm(surf_positions - com, axis=1)
    surf_vectors = (surf_positions - com) / np.array([surf_magnitudes]).T

    #Weight probability of choosing a vector by its length from the center
    surf_probabilities = surf_magnitudes / sum(surf_magnitudes)

    move_vectors_magnitudes_indices = np.random.choice(list(range(len(surf_indices))),
                                                       replace=True,
                                                       size=max_natoms, 
                                                       p=surf_probabilities)

    # Move the atoms
    max_natoms = max(max_natoms, 1)
    natoms_to_move = random.randint(1, max_natoms)
    atom_indices = list(range(len(individual)))
    random.shuffle(atom_indices)  # Using random.shuffle on the indices guarantees no duplicates
    atom_indices = atom_indices[:natoms_to_move]
    for atom_index, move_index in zip(atom_indices, move_vectors_magnitudes_indices):
        vec, mag = surf_vectors[move_index], surf_magnitudes[move_index]
        atom = individual[atom_index]
        atom.x, atom.y, atom.z = com + vec * random.random() * mag

    return None
