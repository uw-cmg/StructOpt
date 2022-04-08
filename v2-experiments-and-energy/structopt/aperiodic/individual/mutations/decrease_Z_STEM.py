import random

import numpy as np
from scipy.ndimage import filters
from ase.data import atomic_numbers, chemical_symbols

from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

def decrease_Z_STEM(individual, STEM_parameters, filter_size=0.5,
                    move_cutoff=0.5, max_cutoff=0.5):
    """Increaes the proton count of an atom to match a STEM image

    Parameters
    ----------
    STEM_parameters : dict
        Parameters for the STEM calculation. Ideally should be the same as the ones
        used for the STEM fitness/relaxation
    filter_size : float
        Filter size for choosing local maximum in the picture. Filter size is equal
        to average_bond_length * resolution * filter_size.
    move_cutoff : float
        The search radius for selecting an atom to move near a high intensity point.
        Defaults to the average bond distance
    """

    module = STEM({'kwargs': STEM_parameters})
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target
    max_max = np.max(contrast)
    min_min = np.min(contrast)

    ###################################
    ## Code for testing the contrast ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(contrast, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))
    # plt.show()
    # import sys; sys.exit()

    # Find a list of local maximum
    cutoff = get_avg_radii(individual) * 2 * 1.1
    move_cutoff *= cutoff
    resolution = module.parameters['kwargs']['resolution']
    size = cutoff * resolution * filter_size

    data_max = filters.maximum_filter(contrast, size=size)
    maxima = ((contrast == data_max) & (contrast > max_max * max_cutoff))
    if len(maxima) == 0:
        return False

    max_coords = np.argwhere(maxima)
    max_xys = (max_coords[:,::-1] - [x_shift, y_shift]) / resolution
    max_intensities = np.asarray([data_max[tuple(coord)] for coord in max_coords])
    max_intensities = np.absolute(max_intensities)
    max_intensities /= sum(max_intensities)

    ###################################
    ## Code for testing the max find ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots(num=1)
    # fig.colorbar(ax.pcolormesh(maxima, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))

    # fig, ax = plt.subplots(num=2)
    # fig.colorbar(ax.pcolormesh(data_max, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))
    
    # plt.show()
    # print(len(max_intensities))
    # import sys; sys.exit()

    # Get atoms associated with each max and max column
    xys = individual.get_positions()[:, :2]

    max_dists = np.expand_dims(xys, 0) - np.transpose(np.expand_dims(max_xys, 0), (1, 0, 2))
    max_dists = np.linalg.norm(max_dists, axis=2)
    max_column_indices = [np.where(dists < move_cutoff)[0] for dists in max_dists]

    if np.size(max_column_indices) == 0:
        return False

    # Elimaxate columns that cannot be "improved" by a permutation
    syms = np.asarray(individual.get_chemical_symbols())
    unique_syms = np.unique(syms)
    unique_nums = [atomic_numbers[sym] for sym in unique_syms]
    max_sym = unique_syms[np.argmax(unique_nums)]
    min_sym = unique_syms[np.argmin(unique_nums)]

    for i, indices in reversed(list(enumerate(max_column_indices))):
        if all(syms[indices] == min_sym):
            max_column_indices = np.delete(max_column_indices, i, axis=0)
            max_intensities = np.delete(max_intensities, i)
            max_intensities /= np.sum(max_intensities)

    # Pick a max column and max column based on their intensities 
    if np.size(max_column_indices) == 0:
        return False
    max_column_indices = max_column_indices[np.random.choice(np.arange(len(max_intensities)), p=max_intensities)]

    # Take out elements with atomic numbers that cannot be increased
    max_column_numbers = [atomic_numbers[sym] for sym in syms[max_column_indices]]
    max_column_indices = [i for i, Z in zip(max_column_indices, max_column_numbers) if Z > np.min(unique_nums)]
    max_column_numbers = [atomic_numbers[sym] for sym in syms[max_column_indices]]

    # Choose a random indice
    max_index = np.random.choice(max_column_indices)
    max_number = atomic_numbers[syms[max_index]]
    new_symbol = chemical_symbols[np.random.choice([n for n in unique_nums if n < max_number])]

    # Switch the atomic symbols
    individual[max_index].symbol = new_symbol

    return
