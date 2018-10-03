from ase.cluster import Octahedron, Icosahedron, Decahedron
from ase.visualize import view
from ase.io import write

atoms = Octahedron('Au', length=5, cutoff=2)
com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Au55-cuboctahedron.xyz', atoms)

atoms = Decahedron('Au', p=3, q=3, r=0)
com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Au55-decahedron.xyz', atoms)

atoms = Icosahedron('Au', noshells=3)
com = atoms.get_center_of_mass()
atoms.set_cell([20, 20, 20])
atoms.center()

write('Au55-icosahedron.xyz', atoms)
