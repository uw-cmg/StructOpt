import random
import numpy as np

from ase import Atom, Atoms
from ase.visualize import view

from structopt.common.crossmodule import NeighborList
from structopt.tools import random_three_vector

def twist(individual, max_radius=0.90):
    """Splits the particle randomly in half and rotates one half.
    
    Parameters
    ----------
    individual : structopt.Individual object
        Individual to be mutated
    max_natoms : float
        That maximum relative distance from the center of the particle 
        the twist is initiated
    """

    if not len(individual):
        return None

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

    # Choose a random point inside the particle
    i = np.random.choice(list(range(len(surf_vectors))), p=surf_probabilities)
    center = com + surf_vectors[i] * surf_magnitudes[i] * random.random() * max_radius

    atoms = individual.copy()
    atoms.translate(-center)
    v = random_three_vector()
    a = random.uniform(30, 180) * np.pi / 180
    atoms.rotate(v=v, a=a, center=(0, 0, 0))
    new_pos_indices = []
    top_atoms = Atoms()
    for i, atom in enumerate(atoms):
        if atom.z > 0:
            new_pos_indices.append(i)
            top_atoms.extend(atom)

    xs, ys, zs = top_atoms.get_positions().T
    x = (max(xs) + min(xs)) * 0.5
    y = (max(ys) + min(ys)) * 0.5
    
    top_atoms.rotate('z', random.uniform(0, np.pi), center=(x, y, 0))
    top_atoms.rotate(v=v, a=-a, center=(0, 0, 0))
    top_atoms.translate(center)
    rotated_positions = top_atoms.get_positions()
    positions = individual.get_positions()
    for i, index in enumerate(new_pos_indices):
        positions[index] = rotated_positions[i]
                     
    individual.set_positions(positions)

    return None
