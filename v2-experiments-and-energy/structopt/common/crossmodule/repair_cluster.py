import random
import numpy as np

from ase import Atoms

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii

def repair_cluster(individual, target_atomlist, surf_CN=11):
    '''Function for repairing a cluster to the correct
    stoichiometry after a crossover'''

    syms = individual.get_chemical_symbols()
    target_atomlist = dict(target_atomlist)
    n_tot = sum([target_atomlist[sym] for sym in target_atomlist])
    atomlist = {sym: syms.count(sym) for sym in set(syms)}

    # If some atoms are missing, make sure they are in the atomlist
    for atom in target_atomlist:
        if atom not in atomlist:
            atomlist[atom] = 0

    # The function depends on how many atoms we have
    n_diff = len(individual) - n_tot
    if n_diff > 0:
        delete_atoms(individual, n_diff, surf_CN)
    elif n_diff < 0:
        add_atoms(individual, n_diff, atomlist, target_atomlist, surf_CN)

    repair_stoichiometry(individual, target_atomlist)
        
def delete_atoms(individual, n_diff, surf_CN):
    '''Function that deletes atoms to make the total atoms correct'''

    syms = individual.get_chemical_symbols()
    
    CNs = CoordinationNumbers(individual)
    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN]
    CNs = np.array(CNs)[surf_indices]
    count_CNs = np.bincount(CNs)
    surf_prob_factors = [1.0 / count_CNs[CN] for CN in CNs]
    surf_probs = 2.0 ** (surf_CN + 1 - CNs)
    surf_probs /= surf_prob_factors
    surf_probs /= np.sum(surf_probs)
    n_delete = np.random.choice(surf_indices, n_diff, False, surf_probs)
    for i in reversed(sorted(n_delete)):
        individual.pop(i)

    return

def add_atoms(individual, n_diff, atomlist, target_atomlist, surf_CN):
    '''Function that adds atoms to make the total atoms correct'''

    # Get a list of symbols for adding
    syms = individual.get_chemical_symbols()
    
    difflist = {sym: atomlist[sym] - target_atomlist[sym] for sym in atomlist}
    add_syms = []
    for sym in difflist:
        if difflist[sym] < 0:
            add_syms += [sym] * -difflist[sym]

    random.shuffle(add_syms)
    add_syms = add_syms[:-n_diff]

    # Get a list of surface atoms to add atoms to 
    CNs = CoordinationNumbers(individual)
    surf_indices = [i for i, CN in enumerate(CNs) if CN <= surf_CN]
    CNs = np.array(CNs)[surf_indices]
    count_CNs = np.bincount(CNs)
    surf_prob_factors = [1.0 / count_CNs[CN] for CN in CNs]
    surf_probs = 2.0 ** CNs
    surf_probs /= surf_prob_factors
    surf_probs /= np.sum(surf_probs)
    if len(individual) < -n_diff:
        add_indices = np.random.choice(surf_indices, -n_diff, True, surf_probs)
    else:
        add_indices = np.random.choice(surf_indices, -n_diff, False, surf_probs)

    # Get positions to add atoms to 
    surf_positions = individual.get_positions()[add_indices]
    COP = individual.get_positions().mean(axis=0)
    vec = surf_positions - COP
    vec /= np.array([np.linalg.norm(vec, axis=1)]).T
    cutoff = get_avg_radii(individual) * 2
    surf_positions = surf_positions + vec * cutoff

    individual.extend(Atoms(symbols=add_syms, positions=surf_positions))

    return

def repair_stoichiometry(individual, target_atomlist):
    syms = individual.get_chemical_symbols()
    atomlist = {sym: syms.count(sym) for sym in set(syms)}

    # If some atoms are missing, make sure they are in the atomlist
    for atom in target_atomlist:
        if atom not in atomlist:
            atomlist[atom] = 0

    if atomlist == target_atomlist:
        return
    
    difflist = {sym: atomlist[sym] - target_atomlist[sym] for sym in atomlist}
    add_dict = {sym: np.absolute(difflist[sym]) for sym in difflist if difflist[sym] < 0}
    del_dict = {sym: np.absolute(difflist[sym]) for sym in difflist if difflist[sym] > 0}

    add_indices = []
    for sym in del_dict:
        new_indices = [i for i, s in enumerate(individual.get_chemical_symbols()) if s == sym]
        random.shuffle(new_indices)
        add_indices += new_indices[:del_dict[sym]]

    add_symbols = []
    for sym in add_dict:
        add_symbols += [sym] * add_dict[sym]

    add_indices = np.sort(add_indices)
    random.shuffle(add_symbols)
    syms = np.asarray(syms)
    syms[add_indices] = add_symbols
    individual.set_chemical_symbols(syms)

    return
