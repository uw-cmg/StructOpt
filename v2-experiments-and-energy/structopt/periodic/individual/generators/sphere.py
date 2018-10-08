import json
import random
import numpy as np
from ase import Atom, Atoms
from ase.visualize import view
from ase.data import atomic_numbers, reference_states

from structopt.tools import random_three_vector
from structopt.common.crossmodule import get_particle_radius

def sphere(atomlist, fill_factor=0.74, radius=None, cell=None):
    """Generates a random sphere of particles given an
    atomlist and radius. If radius is None, one is 
    automatically estimated. min_dist and tries_b4_expand
    are parameters that govern how stricly the proximity
    of atoms are enforced.
    """

    if radius is None:
        radius = get_particle_radius(atomlist, fill_factor)

    # Create a list of random order of the atoms
    chemical_symbols = []
    for atom in atomlist:
        chemical_symbols += [atom[0]] * atom[1]

    random.shuffle(chemical_symbols)

    unit_vec = np.array([random_three_vector() for i in range(len(chemical_symbols))])
    D = radius * np.random.sample(size=len(chemical_symbols)) ** (1.0/3.0)
    positions = np.array([D]).T * unit_vec

    indiv = Atoms(symbols=chemical_symbols, positions=positions)

    if cell is not None:
        indiv.set_cell(cell)
        cell_center = np.sum(indiv.get_cell(), axis=1) / 2.0
        indiv.translate(indiv.get_center_of_mass() + cell_center)
        indiv.set_pbc(True)

    return indiv
