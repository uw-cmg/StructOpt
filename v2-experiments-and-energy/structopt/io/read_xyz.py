import ase.io


def read_xyz(filename, index=None, format=None, **kwargs):
    """Reads an xyz file into an ASE Atoms object and returns it."""
    atoms = ase.io.read(filename, index, format, **kwargs)
    f = open(filename, 'r')
    f.readline()  # natoms
    line = f.readline()  # comment; contains unit cell?
    try:
        line = line.strip().split()[:3]
        x, y, z = float(line[0]), float(line[1]), float(line[2])
        atoms.set_cell([[x, 0, 0], [0, y, 0], [0, 0, z]])
    except Exception as error:
        print(error)
    atoms.set_pbc(True)
    atoms.wrap()
    atoms.center()
    return atoms

