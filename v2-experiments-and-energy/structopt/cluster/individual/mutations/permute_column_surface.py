import random

import numpy as np
from scipy.ndimage import filters, sobel

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

from ase.io import write

def permute_column_surface(individual, STEM_parameters, filter_size=0.5,
                           column_cutoff=0.5):
    """Permutes a column by shifting atoms up and down and filling defects.
    """

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
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(image_max, cmap=cm.viridis))
    # ax.set_aspect('equal', 'box')
    # plt.show()
    # import sys; sys.exit()

    # Pick a random column and find atoms near the column
    i = random.randint(0, len(column_xys) - 1)
    column_xy = column_xys[i]

    xys = individual.get_positions()[:,:2]
    dists = np.linalg.norm(column_xy - xys, axis=1)
    indices = np.arange(len(individual))[dists < column_cutoff]

    if len(indices) < 2:
        return False

    z_positions = individual.get_positions()[indices, 2]
    top_atom = indices[np.argmax(z_positions)]
    bot_atom = indices[np.argmin(z_positions)]

    # Move the top atom the bottom atom or vice versa
    if random.random() < 0.5:
        individual[top_atom].z = individual[bot_atom].z - avg_bond_length
    else:
        individual[bot_atom].z = individual[top_atom].z + avg_bond_length

    return
