import random
import numpy as np

from structopt.common.crossmodule import NeighborList

def move_atoms_group(individual, max_natoms=0.20):
    """Randomly moves atoms within a cluster.
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be moved is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be moved is max_natoms
                                   default: 0.20
    """

    if not len(individual):
        return None

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

    # Weight probability of choosing a vector by its length from the center
    surf_probabilities = surf_magnitudes / sum(surf_magnitudes)

    # Choose a random point inside the particle and find nearest neighbors to that point
    i = np.random.choice(list(range(len(surf_vectors))), p=surf_probabilities)
    center = com + surf_vectors[i] * surf_magnitudes[i] * random.random()
    dists = np.linalg.norm(positions - center, axis=1)
    indices_dists = [[i, dist] for i, dist in enumerate(dists)]
    indices_dists.sort(key=lambda i: i[1])
    max_natoms = max(max_natoms, 1)
    natoms_to_move = random.randint(1, max_natoms)

    # Choose another random point inside the particle and move cluster to that point
    i = np.random.choice(list(range(len(surf_vectors))), p=surf_probabilities)
    move = com + surf_vectors[i] * surf_magnitudes[i] * random.random() - center
    for index, dist in indices_dists[:natoms_to_move]:
        positions[index] += move

    individual.set_positions(positions)

    return None
