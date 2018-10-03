import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, reference_states

def get_avg_radii(atomlist):
    """Returns the average atomic radius of a list of
    atoms. The radius is the radius of the close packed sphere
    in a given crystal structure.

    Parameters
    ----------
    atomlist: list or atoms object
        An N x M list where N is the number of unique atoms and M are
        the attributes of the atom. The first and second index of each
        N list must be the chemical symbol and number, respectively.

    Returns
    -------
    out: float
         The average atomic radii of atomlist    
    """

    # Get the average atomic radii of close packed atoms
    if isinstance(atomlist, Atoms):
        chemical_symbols = atomlist.get_chemical_symbols()
        unique_symbols = set(chemical_symbols)
        atomlist = [[symbol, chemical_symbols.count(symbol)]
                    for symbol in unique_symbols]

    n_tot = sum([atom[1] for atom in atomlist])
    r = 0
    for atom in atomlist:
        n = atom[1]
        conc = float(n)/float(n_tot)
        atomic_number = atomic_numbers[atom[0]]
        struct = reference_states[atomic_number]['symmetry']
        if struct == 'fcc':
            a = reference_states[atomic_number]['a']
            r += conc * np.linalg.norm([a, a]) / 4.0
        elif struct == 'bcc':
            a = reference_states[atomic_number]['a']
            r += conc * np.linalg.norm([a, a, a]) / 4.0
        elif struct == 'hcp':
            a = reference_states[atomic_number]['a']
            r += conc * a
        else:
            raise NotImplementedError('{} structure not supported yet'.format(struct))

    return r
