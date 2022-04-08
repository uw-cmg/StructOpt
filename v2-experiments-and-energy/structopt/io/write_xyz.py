def write_xyz(fileobj, atoms, comment='', append=False):
    """Writes xyz files from an Individual object. Adapted from ase.io.xyz."""

    if isinstance(fileobj, str) and append:
        fileobj = open(fileobj, 'a')
    else:
        fileobj = open(fileobj, 'w')

    symbols = atoms.get_chemical_symbols()
    natoms = len(symbols)

    fileobj.write('{}\n'.format(natoms))
    fileobj.write('{}\n'.format(str(comment)))
    
    atom_positions = atoms.get_positions().copy()
    box_size = atoms.cell[0][0]
    for i, coor in enumerate(atom_positions):
        for j, num in enumerate(coor):
            if num < 0:
                coor[j] = num+box_size
            elif num > box_size:
                coor[j] = num - box_size
            atom_positions[i]=coor
    for s, (x, y, z) in zip(symbols, atoms.get_positions()):
        fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))

