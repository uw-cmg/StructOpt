import random
import numpy as np


def swap_species(individual, max_natoms=0.20):
    """Randomly swaps the species of atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be swapped is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be swapped is max_natoms
                                   if the number of atoms to be swapped is (or evaluates to) an odd integer, it is rounded down to an even integer
                                   max_natoms corresponds to the maximum number of atoms whose species will change
                                   default: 0.20
    """
    if len(individual):
        cell_max = np.maximum.reduce(individual.get_cell())
        cell_min = np.minimum.reduce(individual.get_cell())

        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        max_natoms = max(max_natoms, 2)
        natoms_to_move = random.randint(2, max_natoms)
        if natoms_to_move % 2:  # Make even
            natoms_to_move -= 1
        atom_indices = list(range(len(individual)))
        random.shuffle(atom_indices)  # Using random.shuffle on the indices guarantees no duplicates
        atom_indices = atom_indices[:natoms_to_move]
        # Convert e.g. [1, 2, 3, 4, 5, 6, 7, 8] to [(1, 2), (3, 4), (5, 6), (7, 8)]
        pairs = [(atom_indices[2*i], atom_indices[2*i+1]) for i in range(natoms_to_move//2)]
        for index1, index2 in pairs:
            atom1 = individual[index1]
            atom2 = individual[index2]
            atom1.symbol, atom2.symbol = atom2.symbol, atom1.symbol
    return None

