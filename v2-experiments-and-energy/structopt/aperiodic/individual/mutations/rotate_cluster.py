import random
import numpy as np
from ase import Atom, Atoms

from structopt.common.crossmodule import NeighborList

def rotate_cluster(individual, max_natoms=0.20):
    """Chooses a random number of atoms nearest to a random point in
    the cluster. These atoms are then rotated randomly around this point

    Parameters
    ----------
    individual : Individual
        An individual object
    max_natoms : float
        The fraction of the total atoms to rotate
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
    point = com + surf_vectors[i] * surf_magnitudes[i] * random.random()
    atom = Atom('Si', point)
    nearest_indices = individual.get_nearest_atom_indices(atom_index=atom.index, count=max_natoms)

    # Extract out the atoms to be rotated. Popping changes the indices
    # so pop them out highest to lowest indice
    nearest_indices = np.sort(nearest_indices)[::-1]
    atoms = Atoms()
    for ind in nearest_indices:
        atoms.append(individual.pop(ind))

    axis = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
    angle = random.uniform(30, 180)
    atoms.rotate(axis, a=angle, center='COM', rotate_cell=False)
    individual.extend(atoms)

    return None
