import random

import numpy as np
from scipy.ndimage import filters, sobel

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

from ase.io import write

def enrich_bulk_column(individual, STEM_parameters, filter_size=0.5,
                       column_cutoff=0.5, species=None, surf_CN=11):
    """This mutation randomly does an atom swap within a column of atoms"""

    CNs = CoordinationNumbers(individual)
    
    module = STEM({'kwargs': STEM_parameters})
    module.generate_target()

    image = module.get_image(individual)

    # Find the xy coordinates of the columns in the image
    avg_bond_length = get_avg_radii(individual) * 2
    cutoff = avg_bond_length * 1.1
    column_cutoff *= cutoff
    resolution = module.parameters['kwargs']['resolution']        
    size = cutoff * resolution * filter_size

    image_max = filters.maximum_filter(image, size=size)
    columns = ((image == image_max) & (image > 0.01))
    column_coords = np.argwhere(columns)
    column_xys = column_coords[:, ::-1] / resolution

    if len(column_xys) < 2:
        return False

    ##########################################################
    ### This code is for checking the local maximum finder ###
    ##########################################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots(num=1)
    # fig.colorbar(ax.pcolormesh(image_max, cmap=cm.viridis))
    # ax.set_aspect('equal', 'box')

    # fig, ax = plt.subplots(num=2)
    # fig.colorbar(ax.pcolormesh(columns, cmap=cm.viridis))
    # ax.set_aspect('equal', 'box')
    
    # plt.show()
    # import sys; sys.exit()

    # Get the symbols in each column with > 1 type of atom and a non-species site
    # at the surface available to be switched
    xys = individual.get_positions()[:, :2]
    dists_to_columns = np.expand_dims(xys, 0) - np.transpose(np.expand_dims(column_xys, 0), (1, 0, 2))
    dists_to_columns = np.linalg.norm(dists_to_columns, axis=2)
    all_column_indices = [np.where(dists < column_cutoff)[0] for dists in dists_to_columns]
    syms = individual.get_chemical_symbols()
    if species is None:
        unique_syms = np.unique(syms)
        counts = [syms.count(sym) for sym in unique_syms]
        species = unique_syms[np.argmin(counts)]
    syms = np.asarray(syms)

    column_indices = []
    for indices in all_column_indices:
        column_syms = syms[indices]
        unique_syms = np.unique(column_syms)
        column_CNs = CNs[indices]
        species_CNs = [CN for sym, CN in zip(column_syms, column_CNs) if sym == species]
        other_CNs = [CN for sym, CN in zip(column_syms, column_CNs) if sym != species]
        diff_CNs = np.expand_dims(species_CNs, 0) - np.expand_dims(other_CNs, 0).T
        if len(unique_syms) > 1 and not (diff_CNs >= 0).all():
            column_indices.append(indices)

    if len(column_indices) == 0:
        return False
    
    # Pick a random column and calculate coordination numbers
    column_indices = column_indices[random.randint(0, len(column_indices) - 1)]
    column_CNs = CNs[column_indices]
    species_indices_CNs = [[index, CN] for index, CN in zip(column_indices, column_CNs) if syms[index] == species]
    other_indices_CNs = [[index, CN] for index, CN in zip(column_indices, column_CNs) if syms[index] != species]

    species_indices, species_CNs = zip(*species_indices_CNs)
    other_indices, other_CNs = zip(*other_indices_CNs)

    # Probabilistically select moves that decrease the CN of the enriched species
    diff_CNs = np.expand_dims(species_CNs, 0) - np.expand_dims(other_CNs, 0).T
    moves = np.argwhere(diff_CNs < 0)
    diffs = [diff_CNs[i, j] for i, j in moves]
    diff_counts = {diff: diffs.count(diff) for diff in set(diffs)}
    moves_p = 2.0 ** (np.max(diffs) - np.array(diffs) + 1)
    moves_p /= np.array([diff_counts[diff] for diff in diffs])
    moves_p /= np.sum(moves_p)
    move = moves[np.random.choice(range(len(moves)), p=moves_p)]
    species_index, other_index = species_indices[move[1]], other_indices[move[0]]
    other_symbol = syms[other_index]

    # Apply the mutation
    individual[species_index].symbol = other_symbol
    individual[other_index].symbol = species

    return
