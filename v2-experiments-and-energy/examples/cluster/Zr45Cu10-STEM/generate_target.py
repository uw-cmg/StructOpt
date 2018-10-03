import random

from ase.cluster import Octahedron, Icosahedron, Decahedron
from ase.visualize import view
from ase.io import write

a = 3.23 * 2 ** 0.5

atoms = Octahedron('Zr', latticeconstant=a, length=5, cutoff=2)

indices = list(range(len(atoms)))
random.shuffle(indices)
for i in indices[:10]:
    atoms[i].symbol = 'Cu'

com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Zr45Cu10-cuboctahedron.xyz', atoms)

atoms = Decahedron('Zr', latticeconstant=a, p=3, q=3, r=0)

indices = list(range(len(atoms)))
random.shuffle(indices)
for i in indices[:10]:
    atoms[i].symbol = 'Cu'

com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Zr45Cu10-decahedron.xyz', atoms)

atoms = Icosahedron('Zr', latticeconstant=a, noshells=3)

indices = list(range(len(atoms)))
random.shuffle(indices)
for i in indices[:10]:
    atoms[i].symbol = 'Cu'

com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Zr45Cu10-icosahedron.xyz', atoms)

view(atoms)
