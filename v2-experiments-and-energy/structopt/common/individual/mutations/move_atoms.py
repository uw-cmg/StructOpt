import random
import numpy as np


def move_atoms(individual, max_natoms=0.10):
    """Randomly moves atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be moved is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be moved is max_natoms
                                   default: 0.20
    """
    if len(individual):
        cell_max = np.maximum.reduce(individual.get_cell())
        cell_min = np.minimum.reduce(individual.get_cell())

        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        max_natoms = max(max_natoms, 1)
        natoms_to_move = random.randint(1, max_natoms)
        atom_indices = list(range(len(individual)))
        random.shuffle(atom_indices)  # Using random.shuffle on the indices guarantees no duplicates
        atom_indices = atom_indices[:natoms_to_move]
        for index in atom_indices:
            atom = individual[index]
            atom.x, atom.y, atom.z = (random.uniform(cell_min[0], cell_max[0]),
                                      random.uniform(cell_min[1], cell_max[1]),
                                      random.uniform(cell_min[2], cell_max[2]))
    return None

