import random
import numpy as np
from ase import Atom, Atoms


def rotate_cluster(individual, max_natoms=0.20):
    """Randomly rotates a random cluster of atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be rotated is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be rotated is max_natoms
                                   default: 0.20
    """
    if len(individual):
        cell_max = np.maximum.reduce(individual.get_cell())
        cell_min = np.minimum.reduce(individual.get_cell())

        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        max_natoms = max(max_natoms, 1)
        natoms_to_rotate = random.randint(1, max_natoms)

        point = (random.uniform(cell_min, cell_max),
                 random.uniform(cell_min, cell_max),
                 random.uniform(cell_min, cell_max))

        atom = Atom('Si', point)
        nearest_indices = individual.get_nearest_atom_indices(atom_index=atom.index, count=natoms_to_rotate)

        # Extract out the atoms to be rotated. Popping changes the indices
        # so pop them out highest to lowest indice
        nearest_indices = np.sort(nearest_indices)[::-1]
        atoms = Atoms()
        for ind in nearest_indices:
            atoms.append(individual.pop(ind))

        axis = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
        angle = random.uniform(30, 180)
        atoms.rotate(a=angle, v=axis, center='COM', rotate_cell=False)
        individual.extend(atoms)
    return None

