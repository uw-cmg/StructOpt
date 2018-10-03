import random
import numpy as np
from ase import Atom, Atoms

def rotate_atoms(individual, max_natoms=0.20):
    """Randomly rotates a number of random atoms within the individual (in place).
    
    Args:
        individual (Individual): an individual
        max_natoms (float or int): if float, the maximum number of atoms that will be rotated is max_natoms*len(individual)
                                   if int, the maximum number of atoms that will be rotated is max_natoms
                                   default: 0.20
    """
    if len(individual):
        if isinstance(max_natoms, float):
            max_natoms = int(len(individual)*max_natoms)
        max_natoms = max(max_natoms, 1)
        natoms_to_rotate = random.randint(1, max_natoms)
        atom_indices = list(range(len(individual)))
        random.shuffle(atom_indices)  # Using random.shuffle on the indices guarantees no duplicates
        atom_indices = atom_indices[:natoms_to_rotate]

        # Extract out the atoms to be rotated
        atom_indices.sort(reverse=True)
        atoms = Atoms()
        for ind in atom_indices:
            atoms.append(individual.pop(ind))        
        
        axis = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
        angle = random.uniform(30, 180)
        atoms.rotate(axis, a=angle, center='COM', rotate_cell=False)
        individual.extend(atoms)
    return None

