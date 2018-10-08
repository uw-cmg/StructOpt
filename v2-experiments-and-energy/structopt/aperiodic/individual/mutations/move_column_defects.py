import random
import numpy as np

from ase.data import chemical_symbols
from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii

def move_column_defects(individual, cutoff=0.2, CN_factor=1.1):
    """Calculates the error per column of atoms in the z-direction"""

    avg_radii = get_avg_radii(individual)
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
    for indices in column_indices:
        zs = pos[indices][:,2]
        top_indices.append(indices[np.argmax(zs)])
        bot_indices.append(indices[np.argmin(zs)])

        zs = np.sort(zs)
        if len(zs) == 1:
            avg_bond_lengths.append(np.nan)
        else:
            avg_bond_length = np.average([zs[i+1] - zs[i] for i in range(len(zs)-1)])
            avg_bond_lengths.append(avg_bond_length)

    avg_bond_lengths = np.array(avg_bond_lengths)
    avg_bond_length = np.average(avg_bond_lengths[np.invert(np.isnan(avg_bond_lengths))])
    avg_bond_lengths[np.isnan(avg_bond_lengths)] = avg_bond_length
    avg_bond_lengths = np.append(np.zeros((len(avg_bond_lengths), 2)), np.expand_dims(avg_bond_lengths, 1),  axis=1)

    # Create a list of new surface sites
    bot_new_pos = pos[np.array(bot_indices)] - avg_bond_lengths * 0.95
    top_new_pos = pos[np.array(top_indices)] + avg_bond_lengths * 0.95
    
    # Calculate CNs of old sites and potential new sites
    CNs = CoordinationNumbers(individual, factor=CN_factor)
    CN_cutoff = avg_radii * 2 * CN_factor
    top_CNs = CNs[np.array(top_indices)]
    bot_CNs = CNs[np.array(bot_indices)]

    bot_new_vecs = np.transpose(np.expand_dims(bot_new_pos, 0), [1, 0, 2]) - np.expand_dims(pos, 0)
    bot_new_dists = np.linalg.norm(bot_new_vecs, axis=2)
    bot_new_bonds = (bot_new_dists < CN_cutoff).astype(int)
    bot_new_CNs = np.sum(bot_new_bonds, axis=1)

    top_new_vecs = np.transpose(np.expand_dims(top_new_pos, 0), [1, 0, 2]) - np.expand_dims(pos, 0)
    top_new_dists = np.linalg.norm(top_new_vecs, axis=2)
    top_new_bonds = (top_new_dists < CN_cutoff).astype(int)
    top_new_CNs = np.sum(top_new_bonds, axis=1)

    # Choose moves based on difference in coordination numbers
    move_indices = top_indices + bot_indices
    move_new_pos = np.append(bot_new_pos, top_new_pos, axis=0)
    top_to_bot_CNs = bot_new_CNs - top_CNs
    bot_to_top_CNs = top_new_CNs - bot_CNs
    move_CNs = list(np.append(top_to_bot_CNs, bot_to_top_CNs))
    move_CN_counts = {CN: move_CNs.count(CN) for CN in set(move_CNs)}
    move_probs = [2.0 ** CN / move_CN_counts[CN] for CN in move_CNs]
    move_probs = np.array(move_probs) / sum(move_probs)

    move_indices_i = np.random.choice(range(len(move_indices)), p=move_probs)
    move_index = move_indices[move_indices_i]
    move_new_pos = move_new_pos[move_indices_i]
    
    individual.positions[move_index] = move_new_pos

    return
