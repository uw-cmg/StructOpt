from ase.calculators.lammpsrun import Prism


def write_data(filename, individual):
    """Function for writing the atom positions in a seperate file"""

    individual.wrap()
    individual.center()
    prism = Prism(individual.get_cell())

    with open(filename, 'w') as f:
        f.write('{} (written by StructOpt) \n\n'.format(f.name))
        symbols = individual.get_chemical_symbols()
        n_atoms = len(symbols)
        f.write('{} \t atoms \n'.format(n_atoms))
        species = sorted(set(symbols))
        n_atom_types = len(species)
        f.write('{}  atom types\n'.format(n_atom_types))

        pbc = individual.get_pbc()
        xhi, yhi, zhi, xy, xz, yz = prism.get_lammps_prism_str()
        xyzhis = [xhi, yhi, zhi]
        for index, axis in enumerate(['x','y','z']):
            if pbc[index]:    
                f.write('0.0 {}  {}lo {}hi\n'.format(xyzhis[index], axis, axis))
            else:
                xlo = min([ individual.get_positions()[id][index] for id in range(len(individual.get_positions())) ])
                xhi = max([ individual.get_positions()[id][index] for id in range(len(individual.get_positions())) ])
                f.write('{} {}  {}lo {}hi\n'.format(xlo, xhi, axis, axis))
        
        if prism.is_skewed():
            f.write('{} {} {}  xy xz yz\n'.format(xy, xz, yz))
        
        f.write('\n\n')

        f.write('Atoms \n\n')
        for i, r in enumerate(map(prism.pos_to_lammps_str, individual.get_positions())):
            s = species.index(symbols[i]) + 1
            line = '{:>6} {:>3} {} {} {}\n'
            f.write(line.format(*(i+1, s)+tuple(r)))

