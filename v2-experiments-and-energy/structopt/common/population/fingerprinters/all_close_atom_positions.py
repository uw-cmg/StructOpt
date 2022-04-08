import numpy as np
import scipy.cluster.vq


def all_close_atom_positions(individual1, individual2, rtol=None, atol=1.0):
    """Identifies whether the individuals have the same atom positions within a given tolderance.
    Args:
        individual1 (structopt.common.Individual): The first individual to be compared.
        individual2 (structopt.common.Individual): The second individual to be compared.
        rtol (float): The relative acceptable tolerance between atom positions (see numpy.allclose).
        atol (float): The absolute acceptable tolerance between atom positions (see numpy.allclose).

    Returns:
        bool: True if all atom coordinations are within numpy.allclose's level of acceptability, else False.

    Note:
        If the following equation is element-wise True, then all_close_atom_positions returns True.

        `absolute(individual1.positions - individual2.positions) <= (atol + rtol * absolute(individual2.positions))`
    """

    args = {}
    if rtol is not None:
        args["rtol"] = rtol
    if atol is not None:
        args["atol"] = atol

    a1 = individual1.positions.copy()
    a2 = individual2.positions.copy()
    indexes, dists = scipy.cluster.vq.vq(a1, a2)
    a2sorted = a2[indexes]

    return np.allclose(a1, a2sorted, **args)

