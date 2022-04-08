import random
from ase import Atoms

from structopt.tools import root, single_core, parallel


@single_core
def rotate(individual1, individual2, conserve_composition=True):
    """Rotates the two individuals around their centers of mass,
    splits them in half at the xy-plane, then splices them together.
    Maintains number of atoms.

    Args:
        individual1 (Individual): The first parent
        individual2 (Individual): The second parent
        conserve_composition (bool): default True. If True, conserves composition.

    Returns:
        Individual: The first child
        Individual: The second child

    The children are returned without indicies.
    """
    # Preserve starting conditions of individual
    ind1c = individual1.copy()
    ind2c = individual2.copy()

    # Translate individuals so COM is at (0, 0, 0)
    com1 = ind1c.get_center_of_mass()
    ind1c.translate(-1 * com1)
    com2 = ind2c.get_center_of_mass()
    ind2c.translate(-1 * com2)

    # Select random axis and random angle and rotate individuals
    for _ in range(0, 10):
        rax = random.choice(['x', '-x', 'y', '-y', 'z', '-z'])
        rang = random.random() * 90
        ind1c.rotate(v=rax, a=rang, center='COM', rotate_cell=False)
        # Search for atoms in individual 1 that are above the xy plane
        above_xy_plane1 = Atoms(cell=individual1.get_cell(), pbc=individual1.get_pbc())
        indices1 = []
        for atom in ind1c:
            if atom.z >= 0:
                above_xy_plane1.append(atom)
                indices1.append(atom.index)
        if len(above_xy_plane1) < 2 or len(above_xy_plane1) > len(ind1c):
            # Try again; unrotate ind1c
            ind1c.rotate(v=rax, a=-1*rang, center='COM', rotate_cell=False)
        else:
            break
    ind2c.rotate(v=rax, a=rang, center='COM', rotate_cell=False)

    # Generate `above_xy_plane2`, with the same concentration as `above_xy_plane1` if needed
    above_xy_plane2 = Atoms(cell=individual2.get_cell(), pbc=individual2.get_pbc())
    indices2 = []
    dellist = []
    if conserve_composition:
        symbols = list(set(above_xy_plane1.get_chemical_symbols()))
        # The below list contains atoms from ind2c, whereas `above_xy_plane1` contains atoms from ind1c
        atoms_by_symbol = {sym: [atm for atm in ind2c if atm.symbol == sym] for sym in symbols}
        for atom in above_xy_plane1:
            if len(atoms_by_symbol[atom.symbol]) > 0:
                # Get the atom from `atoms_by_symbol` that has the same type as `atom` and the largest `z` value
                dist = [atom.z for atom in atoms_by_symbol[atom.symbol]]
                pos = dist.index(max(dist))
                above_xy_plane2.append(atoms_by_symbol[atom.symbol][pos])
                indices2.append(atoms_by_symbol[atom.symbol][pos].index)
                del atoms_by_symbol[atom.symbol][pos]
            else:
                dellist.append(atom.index)
        if dellist:
            dellist.sort(reverse=True)
            for atom_index in dellist:
                del above_xy_plane1[atom_index]
                del indices1[atom_index]
    else:
        for atom in ind2c:
            if atom.z >= 0:
                above_xy_plane2.append(atom)
                indices2.append(atom.index)
        while len(above_xy_plane2) < len(above_xy_plane1)-len(dellist):
            # Too many atoms in above_xy_plane1
            dellist.append(random.choice(above_xy_plane1).index)
        if dellist:
            dellist.sort(reverse=True)
            for atom in dellist:
                del above_xy_plane1[atom]
                del indices1[atom]
        dellist = []
        while len(above_xy_plane1) < len(above_xy_plane2)-len(dellist):
            # Too many atoms in above_xy_plane2
            dellist.append(random.choice(above_xy_plane2).index)
        if dellist:
            dellist.sort(reverse=True)
            for atom in dellist:
                del above_xy_plane2[atom]
                del indices2[atom]

    below_xy_plane1 = Atoms()
    below_xy_plane2 = Atoms()
    below_xy_plane2 = Atoms(cell=individual2.get_cell(), pbc=individual2.get_pbc())
    for atom in ind2c:
        if atom.index not in indices2:
            below_xy_plane2.append(atom)
    below_xy_plane1 = Atoms(cell=individual1.get_cell(), pbc=individual1.get_pbc())
    for atom in ind1c:
        if atom.index not in indices1:
            below_xy_plane1.append(atom)


    child1 = above_xy_plane2.copy()
    child1.extend(below_xy_plane1)
    child2 = above_xy_plane1.copy()
    child2.extend(below_xy_plane2)

    # Need to have at least one atom of each specie in atomlist to prevent LAMMPS from erroring
    if not conserve_composition:
        for i in range(len(individual1)):
            atoms1 = [atom for atom in child1 if atom.symbol == individual1[i]]
            atoms2 = [atom for atom in child2 if atom.symbol == individual1[i]]
            if len(atoms1) == 0 and len(atoms2) == 0:
                random.choice(child1).symbol = individual1[i].symbol
                random.choice(child2).symbol = individual1[i].symbol
            elif len(atoms1) == 0 and len(atoms2) != 0:
                del child1[random.randint(0, len(child1))]
                child1.append(random.choice(atoms2))
            elif len(atoms1) != 0 and len(atoms2) == 0:
                del child2[random.randint(0, len(child2))]
                child2.append(random.choice(atoms1))

    # Unrotate and untranslate the children
    child1.rotate(v=rax, a=-1*rang, center='COM', rotate_cell=False)
    child2.rotate(v=rax, a=-1*rang, center='COM', rotate_cell=False)
    child1.translate(com1)
    child2.translate(com2)

    full_child1 = individual1.copy()
    full_child1.clear()
    full_child1.extend(child1)
    full_child2 = individual2.copy()
    full_child2.clear()
    full_child2.extend(child2)
    return full_child1, full_child2

