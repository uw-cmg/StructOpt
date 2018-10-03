import random
import numpy as np
from structopt.common.crossmodule import NeighborList

def poor2rich(individual):
    '''Used for multi-component systems. Swaps atoms A and B
    so that atom A moves from a region with a low number of A-A
    bonds to a high number of A-A bonds.

    Parameters
    ----------
    individual : Individual
        An individual
    '''

    syms = np.asarray(individual.get_chemical_symbols())
    NN_list = NeighborList(individual)
    NNs = [list(syms[i]) for i in NN_list]
    max_CN = max([len(NN) for NN in NNs])

    # Pick an atom in a poor environment
    ns_to_rich = np.array([max_CN - NNs[i].count(sym) for i, sym in enumerate(syms)])

    unique_ns = np.array(list(set(ns_to_rich)))
    prob_poor = 2.0 ** unique_ns
    prob_poor /= np.sum(prob_poor)
    n_to_rich = np.random.choice(unique_ns, p=prob_poor)
    indices_poor = np.where(ns_to_rich == n_to_rich)[0]
    index_poor = np.random.choice(indices_poor)
    symbol_poor = syms[index_poor]

    # Pick another atom in a rich environment
    indices_other = [i for i, sym in enumerate(syms) if sym != symbol_poor]
    ns = np.array([NNs[i].count(symbol_poor) for i in indices_other])
    unique_ns = np.array(list(set(ns)))
    prob_rich = 2.0 ** unique_ns
    prob_rich /= np.sum(prob_rich)
    n = np.random.choice(unique_ns, p=prob_rich)
    indices_rich = np.where(ns == n)[0]
    index_rich = np.random.choice(indices_rich)
    index_rich = indices_other[index_rich]
    symbol_rich = syms[index_rich]

    individual[index_poor].symbol = symbol_rich
    individual[index_rich].symbol = symbol_poor

    return
