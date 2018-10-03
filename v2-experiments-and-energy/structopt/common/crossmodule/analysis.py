import numpy as np
from structopt.common.crossmodule import get_avg_radii

np.seterr(all='ignore')

def CoordinationNumbers(atoms, cutoff=None, factor=1.1):
    """Calculates the coordination number of all atoms based on
    cutoff radius "cutoff". Does not obey periodic boundary conditions.

    Parameters
    ----------
    atoms : ase.Atoms or structopt.Individual object
        The atoms object to be analyzed
    cutoff : float
        The radius to search for neighbors. If cutoff is not
        specified, returns average bond length from a weighted
        average of experimental bond lengths.
    factor : float
        If cutoff is None, nearest neighbor distance is taken as
        two times the average atomic radius. factor is used to
        expand the cutoff by cutoff * factor to ensure python
        numerical behavior doesn't "lose" atoms.

    Output
    ------
    out : list
        A list of coordination numbers, where out[i] corresponds
        to the coordination number of atoms[i]
    """

    if cutoff is None:
        chemical_symbols = atoms.get_chemical_symbols()
        unique_symbols = set(chemical_symbols)
        atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
        cutoff = get_avg_radii(atomlist) * 2 * factor

    pos = np.array([atoms.get_positions()])
    pos_T = np.transpose(pos, [1, 0, 2])
    vecs = pos - pos_T
    dists = np.linalg.norm(vecs, axis=2)

    bonds = ((dists < cutoff).astype(int) * (dists > 0).astype(int))
    CNs = np.sum(bonds, axis=1)

    return CNs

def NeighborList(atoms, cutoff=None, factor=1.1):
    """Calculates the neighbors of all atoms based on
    cutoff radius "cutoff". Does not obey periodic boundary conditions.

    Parameters
    ----------
    atoms : ase.Atoms or structopt.Individual object
        The atoms object to be analyzed
    cutoff : float
        The radius to search for neighbors. If cutoff is not
        specified, returns average bond length from a weighted
        average of experimental bond lengths.
    factor : float
        If cutoff is None, nearest neighbor distance is taken as
        two times the average atomic radius. factor is used to
        expand the cutoff by cutoff * factor to ensure python
        numerical behavior doesn't "lose" atoms.

    Output
    ------
    out : list
        A list of a neighbors. out[i] returns a list of the neighbors
        of atom i.
    """

    if cutoff is None:
        chemical_symbols = atoms.get_chemical_symbols()
        unique_symbols = set(chemical_symbols)
        atomlist = [[symbol, chemical_symbols.count(symbol)] for symbol in unique_symbols]
        cutoff = get_avg_radii(atomlist) * 2 * factor

    pos = np.array([atoms.get_positions()])
    pos_T = np.transpose(pos, [1, 0, 2])
    vecs = pos - pos_T
    dists = np.linalg.norm(vecs, axis=2)

    bonds = ((dists < cutoff).astype(int) * (dists > 0).astype(int))

    neighbors = np.resize(np.arange(len(atoms)), (len(atoms), len(atoms)))
    neighbors = neighbors / bonds.astype(float)

    neighbors = [a[~(np.isnan(a) | np.isinf(a))].astype(int) for a in neighbors]    

    return np.asarray(neighbors)

def NeighborElements(atoms, cutoff=None, factor=1.1):
    """Gives the neighboring elements of each atom
    cutoff radius "cutoff". Does not obey periodic boundary conditions.

    Parameters
    ----------
    atoms : ase.Atoms or structopt.Individual object
        The atoms object to be analyzed
    cutoff : float
        The radius to search for neighbors. If cutoff is not
        specified, returns average bond length from a weighted
        average of experimental bond lengths.
    factor : float
        If cutoff is None, nearest neighbor distance is taken as
        two times the average atomic radius. factor is used to
        expand the cutoff by cutoff * factor to ensure python
        numerical behavior doesn't "lose" atoms.

    Output
    ------
    out : list
        A list of a neighbors. out[i] returns a list of the neighbors
        of atom i.
    """

    syms = np.asarray(atoms.get_chemical_symbols())
    neighbors = NeighborList(atoms, cutoff=None, factor=1.1)
    neighbors = [list(syms[i]) for i in neighbors]
    return neighbors
