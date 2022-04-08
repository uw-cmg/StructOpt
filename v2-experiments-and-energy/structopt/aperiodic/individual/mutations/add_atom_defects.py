import random
import numpy as np

from ase import Atom, Atoms
from structopt.common.crossmodule import get_avg_radii

def add_atom_defects(individual, add_prob=None, cutoff=0.2, CN_factor=1.1):
    """Calculates the error per column of atoms in the z-direction"""

    avg_radii = get_avg_radii(individual)
    avg_bond_length = avg_radii * 2
    cutoff *= avg_radii * 2

    # Organize atoms into columns
    pos = individual.get_positions()
    xys = np.expand_dims(pos[:, :2], 0)
    dists = np.linalg.norm(xys - np.transpose(xys, (1, 0, 2)), axis=2)

    NNs = np.sort(np.argwhere(dists < cutoff))
    column_indices = []
    atoms_to_be_sorted = list(range(len(individual)))
    while len(atoms_to_be_sorted) > 0:
        i = atoms_to_be_sorted[0]
        same_column_indices = np.unique(NNs[NNs[:,0] == i])
        column_indices.append(same_column_indices)
        for j in reversed(sorted(same_column_indices)):
            i_del = atoms_to_be_sorted.index(j)
            atoms_to_be_sorted.pop(i_del)
            NNs = NNs[NNs[:,0] != j]
            NNs = NNs[NNs[:,1] != j]

    # Make a list of the top and bottom atom of each column as well
    # the average bond length of atoms in the column
    top_indices, bot_indices, avg_bond_lengths = [], [], []
    vac_new_pos = []
    for indices in column_indices:
        zs = pos[indices][:,2]
        top_indices.append(indices[np.argmax(zs)])
        bot_indices.append(indices[np.argmin(zs)])

        zs = np.sort(zs)
        if len(zs) == 1:
            avg_bond_lengths.append(np.nan)
            continue
        else:
            avg_length = np.average([zs[i+1] - zs[i] for i in range(len(zs)-1)
                                          if zs[i+1] - zs[i] < avg_bond_length * 1.7])
            avg_bond_lengths.append(avg_length)

        # Check for vacancies in the column
        for i, z in enumerate(zs[:-1]):
            diff = zs[i+1] - z
            if diff < avg_length * 1.5:
                continue
            z += avg_bond_length
            xy = pos[indices][:,:2].mean(axis=0)
            vac_new_pos.append(np.append(xy, z))

    avg_bond_lengths = np.array(avg_bond_lengths)
    avg_bond_length = np.average(avg_bond_lengths[np.invert(np.isnan(avg_bond_lengths))])
    avg_bond_lengths[np.isnan(avg_bond_lengths)] = avg_bond_length
    avg_bond_lengths = np.append(np.zeros((len(avg_bond_lengths), 2)), np.expand_dims(avg_bond_lengths, 1),  axis=1)

    # Create a list of new surface sites
    bot_new_pos = pos[np.array(bot_indices)] - avg_bond_lengths * 0.95
    top_new_pos = pos[np.array(top_indices)] + avg_bond_lengths * 0.95
    if np.size(vac_new_pos) > 0:
        surf_positions = np.concatenate((bot_new_pos, top_new_pos, vac_new_pos), axis=0)
    else:
        surf_positions = np.concatenate((bot_new_pos, top_new_pos), axis=0)

    # Calculate CNs of old sites and potential new sites
    CN_cutoff = avg_radii * 2 * CN_factor
    add_vecs = np.transpose(np.expand_dims(surf_positions, 0), [1, 0, 2]) - np.expand_dims(pos, 0)
    add_dists = np.linalg.norm(add_vecs, axis=2)
    add_bonds = (add_dists < CN_cutoff).astype(int)
    add_CNs = list(np.sum(add_bonds, axis=1))

    # Choose moves based on difference in coordination numbers
    add_CN_counts = {CN: add_CNs.count(CN) for CN in set(add_CNs)}
    add_probs = [2.0 ** CN / add_CN_counts[CN] for CN in add_CNs]
    add_probs = np.array(add_probs) / sum(add_probs)

    new_position = surf_positions[np.random.choice(range(len(surf_positions)), p=add_probs)]

    # Choose the element to add
    if add_prob is None:
        syms = individual.get_chemical_symbols()
        elements = np.unique(syms)
        n = len(syms)
        p = [syms.count(element) / n for element in elements]
    else:
        elements = [key for key in add_prob]
        p = [add_prob[key] for key in elements]

    element = np.random.choice(elements, p=p)
    individual.extend(Atoms([Atom(element, new_position)]))

    return
