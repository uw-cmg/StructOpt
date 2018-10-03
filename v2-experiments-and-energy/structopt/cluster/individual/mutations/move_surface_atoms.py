import random
import numpy as np

from structopt.common.crossmodule import NeighborList
from structopt.common.crossmodule import get_avg_radii

from ase.io import write

def move_surface_atoms(individual, max_natoms=0.2, move_CN=11, surf_CN=11):
    """Randomly moves atoms at the surface to other surface sites

    Parameters
    ----------
    individual : Individual
        The individual object to be modified in place
    max_natoms : float or int
        if float, the maximum number of atoms that will be moved is 
        max_natoms*len(individual)
        if int, the maximum number of atoms that will be moved is max_natoms
        default: 0.20
    move_CN : int
        The coordination number to determine which atoms can move moved. Any 
        atom with coordination number above move_CN will not be moved
    surf_CN : int
        The coordination number to determine which atoms are considered surface
        atoms. Surface atoms are used to estimating new surface sites
    """

    if len(individual) == 0:
        return False

    # Analyze the individual
    NNs = NeighborList(individual)
    CNs = [len(NN) for NN in NNs]
    
    # Get indices of atoms considered to be moved
    move_indices_CNs = [[i, CN] for i, CN in enumerate(CNs) if CN <= move_CN]
    if len(move_indices_CNs) == 0:
        return False
    move_indices_CNs.sort(key=lambda i: i[1])
    move_indices = list(zip(*move_indices_CNs))[0]

    # Get surface sites to move atoms to
    # First get all surface atoms
    positions = individual.get_positions()
    surf_indices_CNs = [[i, CN] for i, CN in enumerate(CNs)
                        if CN <= surf_CN and CN > 2]
    surf_indices_CNs.sort(key=lambda i: i[1])
    surf_indices = list(zip(*surf_indices_CNs))[0]
    surf_positions = np.array([positions[i] for i in surf_indices])

    # Get the average bond length of the particle
    chemical_symbols = individual.get_chemical_symbols()
    unique_symbols = set(chemical_symbols)
    atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
    avg_bond_length = get_avg_radii(atomlist) * 2

    # Choose sites as projections one bond length away from COM
    COM = surf_positions.mean(axis=0)
    vec = surf_positions - COM
    vec /= np.array([np.linalg.norm(vec, axis=1)]).T
    add_positions = surf_positions + vec * avg_bond_length * 0.5

    # Set positions of a fraction of the surface atoms
    if type(max_natoms) is float:
        max_natoms = int(max_natoms * len(move_indices))
    move_natoms = random.randint(0, max_natoms)
    move_indices = move_indices[:move_natoms]
    add_indices = np.random.choice(len(add_positions), len(move_indices), replace=False)
    for move_index, add_index in zip(move_indices, add_indices):
        positions[move_index] = add_positions[add_index]

    individual.set_positions(positions)

    return

