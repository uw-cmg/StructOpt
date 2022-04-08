# from ase.calculators.lammpsrun import Prism
import numpy as np
import decimal as dec

class Prism:
    def __init__(self, cell, pbc=(True,True,True), digits=10):
        """Create a lammps-style triclinic prism object from a cell

        The main purpose of the prism-object is to create suitable 
        string representations of prism limits and atom positions
        within the prism.
        When creating the object, the digits parameter (default set to 10)
        specify the precission to use.
        lammps is picky about stuff being within semi-open intervals,
        e.g. for atom positions (when using create_atom in the in-file), 
        x must be within [xlo, xhi).
        """
        a, b, c = cell
        an, bn, cn = [np.linalg.norm(v) for v in cell]
        
        alpha = np.arccos(np.dot(b, c)/(bn*cn))
        beta  = np.arccos(np.dot(a, c)/(an*cn))
        gamma = np.arccos(np.dot(a, b)/(an*bn))
        
        xhi = an
        xyp = np.cos(gamma)*bn
        yhi = np.sin(gamma)*bn
        xzp = np.cos(beta)*cn
        yzp = (bn*cn*np.cos(alpha) - xyp*xzp)/yhi
        zhi = np.sqrt(cn**2 - xzp**2 - yzp**2)
    
        # Set precision
        self.car_prec = dec.Decimal('10.0') ** \
            int(np.floor(np.log10(max((xhi,yhi,zhi))))-digits)
        self.dir_prec = dec.Decimal('10.0') ** (-digits)
        self.acc = float(self.car_prec)
        self.eps = np.finfo(xhi).eps

        # For rotating positions from ase to lammps
        Apre = np.array(((xhi, 0,   0),
                         (xyp, yhi, 0),
                         (xzp, yzp, zhi)))
        self.R = np.dot(np.linalg.inv(cell), Apre)

        # Actual lammps cell may be different from what is used to create R
        def fold(vec, pvec, i):
            p = pvec[i]
            x = vec[i] + 0.5*p
            n = (np.mod(x, p) - x)/p
            return [float(self.f2qdec(a)) for a in (vec + n*pvec)]

        Apre[1,:] = fold(Apre[1,:], Apre[0,:], 0)
        Apre[2,:] = fold(Apre[2,:], Apre[1,:], 1)
        Apre[2,:] = fold(Apre[2,:], Apre[0,:], 0)

        self.A = Apre
        self.Ainv = np.linalg.inv(self.A)

        if self.is_skewed() and \
                (not (pbc[0] and pbc[1] and pbc[2])):
            raise RuntimeError('Skewed lammps cells MUST have '
                               'PBC == True in all directions!')

    def f2qdec(self, f):
        return dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_DOWN)

    def f2qs(self, f):
        return str(self.f2qdec(f))

    def f2s(self, f):
        return str(dec.Decimal(repr(f)).quantize(self.car_prec, dec.ROUND_HALF_EVEN))

    def dir2car(self, v):
        "Direct to cartesian coordinates"
        return np.dot(v, self.A)

    def car2dir(self, v):
        "Cartesian to direct coordinates"
        return np.dot(v, self.Ainv)

    def fold_to_str(self,v):
        "Fold a position into the lammps cell (semi open), return a tuple of str" 
        # Two-stage fold, first into box, then into semi-open interval
        # (within the given precission).
        d = [x % (1-self.dir_prec) for x in 
             map(dec.Decimal, map(repr, np.mod(self.car2dir(v) + self.eps, 1.0)))]
        return tuple([self.f2qs(x) for x in 
                      self.dir2car(map(float, d))])
        
    def get_lammps_prism(self):
        A = self.A
        return (A[0,0], A[1,1], A[2,2], A[1,0], A[2,0], A[2,1])

    def get_lammps_prism_str(self):
        "Return a tuple of strings"
        p = self.get_lammps_prism()
        return tuple([self.f2s(x) for x in p])

    def pos_to_lammps_str(self, position):
        "Rotate an ase-cell postion to the lammps cell orientation, return tuple of strs"
        return tuple([self.f2s(x) for x in np.dot(position, self.R)])

    def pos_to_lammps_fold_str(self, position):
        "Rotate and fold an ase-cell postion into the lammps cell, return tuple of strs"
        return self.fold_to_str(np.dot(position, self.R))

    def is_skewed(self):
        acc = self.acc
        prism = self.get_lammps_prism()
        axy, axz, ayz = [np.abs(x) for x in prism[3:]]
        return (axy >= acc) or (axz >= acc) or (ayz >= acc)


def write_data(filename, individual, atom_style = 'atomic'):
    """Function for writing the atom positions in a seperate file"""

    assert atom_style in ['atomic', 'charge'], 'atom_style only support atomic or charge'
    individual.wrap()
    individual.center()
    prism = Prism(individual.get_cell())

    with open(filename, 'w') as f:
        f.write('{} (written by StructOpt) \n\n'.format(f.name))
        symbols = individual.get_chemical_symbols()
        n_atoms = len(symbols)
        f.write('{} \t atoms \n'.format(n_atoms))
        species = sorted(set(symbols), reverse=True)
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
            if atom_style == 'atomic':
                line = '{:>6} {:>3} {} {} {}\n'
                f.write(line.format(*(i+1, s)+tuple(r)))
            elif atom_style == 'charge':
                line = '{:>6} {:>3} {} {} {} {}\n'
                '''
                if symbols[i] == 'O':
                    # for comb
                    #q = -1.62
                    # for reax/c
                    q = -0.8
                elif symbols[i] == 'Ti':
                    # for comb
                    #q = 3.24
                    # for reax/c
                    q = 1.6
                else:
                    q = 0.0
                '''
                q = 0.0
                f.write(line.format(*(i+1, s, q)+tuple(r)))
