import random

import numpy as np
from scipy.ndimage import filters
from ase.data import atomic_numbers

from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

def permutation_STEM(individual, STEM_parameters, filter_size=0.5,
                     move_cutoff=0.5, max_cutoff=0.5, min_cutoff=0.5):
    """Moves surface atoms around based on the difference in the target
    and individual STEM image

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

    # Find a list of local maximum and local minimum in the image
    cutoff = get_avg_radii(individual) * 2 * 1.1
    move_cutoff *= cutoff
    resolution = module.parameters['kwargs']['resolution']
    size = cutoff * resolution * filter_size

    data_max = filters.maximum_filter(contrast, size=size)
    maxima = ((contrast == data_max) & (contrast > max_max * max_cutoff))
    if len(maxima) == 0:
        return False
    max_coords = np.argwhere(maxima)
    max_xys = (max_coords[:,::-1] - np.array([[x_shift, y_shift]])) / resolution
    max_intensities = np.asarray([data_max[tuple(coord)] for coord in max_coords])
    max_intensities /= sum(max_intensities)

    ###################################
    ## Code for testing the max find ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(maxima, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))
    # print(len(max_intensities))

    # fig, ax = plt.subplots(num=2)
    # fig.colorbar(ax.pcolormesh(data_max, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))
    
    # plt.show()
    # print(len(max_intensities))
    # import sys; sys.exit()

    data_min = filters.minimum_filter(contrast, size=size)
    minima = ((contrast == data_min) & (contrast < min_min * min_cutoff))
    if len(minima) == 0:
        return False

    min_coords = np.argwhere(minima)
    min_xys = (min_coords[:,::-1] - [x_shift, y_shift]) / resolution
    min_intensities = np.asarray([data_min[tuple(coord)] for coord in min_coords])
    min_intensities = np.absolute(min_intensities)
    min_intensities /= sum(min_intensities)

    ###################################
    ## Code for testing the min find ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots(num=1)
    # fig.colorbar(ax.pcolormesh(minima, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))

    # fig, ax = plt.subplots(num=2)
    # fig.colorbar(ax.pcolormesh(data_min, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * STEM_parameters['resolution']))
    
    # plt.show()
    # print(len(min_intensities))
    # import sys; sys.exit()

    # Get atoms associated with each max and min column
    xys = individual.get_positions()[:, :2]
    max_dists = np.expand_dims(xys, 0) - np.transpose(np.expand_dims(max_xys, 0), (1, 0, 2))
    max_dists = np.linalg.norm(max_dists, axis=2)
    max_column_indices = [np.where(dists < move_cutoff)[0] for dists in max_dists]

    min_dists = np.expand_dims(xys, 0) - np.transpose(np.expand_dims(min_xys, 0), (1, 0, 2))
    min_dists = np.linalg.norm(min_dists, axis=2)
    min_column_indices = [np.where(dists < move_cutoff)[0] for dists in min_dists]

    if np.size(min_column_indices) == 0 or np.size(max_column_indices) == 0:
        return False

    # Eliminate columns that cannot be "improved" by a permutation
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

    for i, indices in reversed(list(enumerate(min_column_indices))):
        if all(syms[indices] == max_sym):
            min_column_indices = np.delete(min_column_indices, i, axis=0)
            min_intensities = np.delete(min_intensities, i)
            min_intensities /= np.sum(min_intensities)

    # Pick a max column and min column based on their intensities 
    if np.size(min_column_indices) == 0 or np.size(max_column_indices) == 0:
        return False

    max_column_indices = max_column_indices[np.random.choice(np.arange(len(max_intensities)), p=max_intensities)]
    min_column_indices = min_column_indices[np.random.choice(np.arange(len(min_intensities)), p=min_intensities)]

    # Pick a move between the two columns based on differences in atomic numbers
    max_column_numbers = [atomic_numbers[sym] for sym in syms[max_column_indices]]    
    min_column_numbers = [atomic_numbers[sym] for sym in syms[min_column_indices]]

    max_column_numbers = np.expand_dims(max_column_numbers, 0)    
    min_column_numbers = np.expand_dims(min_column_numbers, 0).T

    min_max_pairs = np.argwhere(max_column_numbers - min_column_numbers > 0)
    min_index, max_index = min_max_pairs[random.randint(0, len(min_max_pairs) - 1)]
    min_index, max_index = min_column_indices[min_index], max_column_indices[max_index]
    max_symbol = syms[max_index]
    min_symbol = syms[min_index]

    # Switch the atomic symbols
    individual[max_index].symbol = min_symbol
    individual[min_index].symbol = max_symbol

    return
