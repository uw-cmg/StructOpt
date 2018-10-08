import random
from math import atan
import numpy as np
np.seterr(all='ignore')

from ase import Atoms
from structopt.tools import random_three_vector
from structopt.common.crossmodule import CoordinationNumbers

from .fcc import fcc

def fcc_110_twin(atomlist, cell, a, shape=[1, 1, 1], roundness=0.5, alpha=10,
                 rotation=0, center=0.5):

    # Make a bigger particle to ensure child has more than required
    n = int(sum(list(zip(*atomlist))[1]) * 1.1)

    atoms1 = fcc([['Pt', n]], cell, a, shape=shape, orientation='110', angle=rotation)
    rotation2 = rotation + (np.pi - atan(2 ** 0.5) * 2)
    atoms2 = fcc([['Pt', n]], cell, a, shape=shape, orientation='110', angle=rotation2)

    # Now rotate both atoms so the twin plane is parallel to the x axis
    align = rotation - atan(2 ** 0.5)
    atoms1.rotate('z', -align, center='COP')
    atoms2.rotate('z', -align, center='COP')

    # Make a new atoms object from the two atoms.
    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()
    mid_x1 = (np.max(pos1[:,0]) + np.min(pos1[:,0])) * 0.5
    mid_y1 = (np.max(pos1[:,1]) + np.min(pos1[:,1])) * center
    mid_z1 = (np.max(pos1[:,2]) + np.min(pos1[:,2])) * 0.5
    center1 = np.array([mid_x1, mid_y1, mid_z1])

    mid_x2 = (np.max(pos2[:,0]) + np.min(pos2[:,0])) * 0.5
    mid_y2 = (np.max(pos2[:,1]) + np.min(pos2[:,1])) * center
    mid_z2 = (np.max(pos2[:,2]) + np.min(pos2[:,2])) * 0.5
    center2 = np.array([mid_z2, mid_y2, mid_z2])

    dists1 = np.linalg.norm(pos1 - center1, axis=1)
    center1 = pos1[np.argmin(dists1)]
    atoms1.translate(-center1)

    dists2 = np.linalg.norm(pos2 - center2, axis=1)
    center2 = pos2[np.argmin(dists2)]
    atoms2.translate(-center2)

    child = atoms1[atoms1.get_positions()[:,1] <= 0.1] + atoms2[atoms2.get_positions()[:,1] > 0.1]
    child.center()
    child.rotate('z', align, center='COP')

    # Delete any extra atoms on the surface
    cutoff = np.linalg.norm([a, a]) / 4.0 * 2 * 1.1
    CNs = CoordinationNumbers(child, cutoff)
    indices_CNs = list(zip(range(len(CNs)), CNs))
    indices_CNs = sorted(indices_CNs, key=lambda i: i[1])
    assert len(CNs) - sum(list(zip(*atomlist))[1])
    to_del = list(zip(*indices_CNs))[0][:len(CNs) - sum(list(zip(*atomlist))[1])]
    to_del = sorted(to_del)
    for i in reversed(to_del):
        child.pop(i)

    # Make the composition correct
    symbols = []
    for sym, n in atomlist:
        symbols += [sym] * n
    random.shuffle(symbols)
    child.set_chemical_symbols(symbols)

    return child
