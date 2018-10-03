import random

import numpy as np
from scipy.ndimage import filters, sobel

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

from ase.io import write

def permute_column_bulk(individual, STEM_parameters, filter_size=0.5,
                        column_cutoff=0.5):
    """This mutation randomly does an atom swap within a column of atoms"""

    module = STEM(STEM_parameters)
    module.generate_target()

    image = module.get_image(individual)

    # Find the xy coordinates of the columns in the image
    avg_bond_length = get_avg_radii(individual) * 2
    cutoff = avg_bond_length * 1.1
    column_cutoff *= cutoff
    resolution = module.parameters['resolution']        
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

    # Get the symbols in each column with > 1 type of atom
    xys = individual.get_positions()[:, :2]
    dists_to_columns = np.expand_dims(xys, 0) - np.transpose(np.expand_dims(column_xys, 0), (1, 0, 2))
    dists_to_columns = np.linalg.norm(dists_to_columns, axis=2)
    column_indices = [np.where(dists < column_cutoff)[0] for dists in dists_to_columns]
    syms = np.asarray(individual.get_chemical_symbols())
    column_indices = [indices for indices in column_indices if len(np.unique((syms[indices]))) > 1]

    # Shuffle the indices in the atomic column
    column_indices = column_indices[random.randint(0, len(column_indices) - 1)]
    column_symbols = syms[column_indices]
    np.random.shuffle(column_symbols)
    for index, symbol in zip(column_indices, column_symbols):
        individual[index].symbol = symbol
    
    return
