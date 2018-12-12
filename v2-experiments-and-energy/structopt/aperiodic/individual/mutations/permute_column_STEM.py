import random

import numpy as np
from scipy.ndimage import filters, sobel

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

from ase.io import write

def permute_column_STEM(individual, STEM_parameters, filter_size=1,
                        column_cutoff=0.5, max_cutoff=0.5):
    """Permutes a column based on the misalignment of columns between
    the individual and target.

    """

    module = STEM({'kwargs': STEM_parameters})
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target
    max_max = np.max(contrast)
    min_min = np.min(contrast)

    # Find a list of local maximum and local minimum in the image
    cutoff = get_avg_radii(individual) * 2 * 1.1
    column_cutoff *= cutoff
    resolution = module.parameters['kwargs']['resolution']        
    size = cutoff * resolution * filter_size

    image_max = filters.maximum_filter(image, size=size)
    columns = ((image == image_max) & (image > 0.01))
    columns = np.argwhere(columns)
    columns = (columns[:, ::-1] - x_shift, y_shift)    

    data_max = filters.maximum_filter(contrast, size=size)
    maxima = ((contrast == data_max) & (contrast > 0.01))
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
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))
    # plt.show()
    # print(len(max_intensities))
    # import sys; sys.exit()

    data_min = filters.minimum_filter(contrast, size=size)
    minima = ((contrast == data_min) & (contrast < -0.01))
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
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(data_min, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))
    # plt.show()
    # print(len(min_intensities))
    # import sys; sys.exit()

    # Find areas where where max and min are next to each other
    mins = np.array([min_xys])
    maxs = np.transpose(np.array([max_xys]), [1, 0, 2])
    vecs = mins - maxs
    dists = np.linalg.norm(vecs, axis=2)

    bonds = (dists < column_cutoff).astype(int)
    CNs = np.sum(bonds, axis=1)
    max_xys = max_xys[CNs > 0]
    max_intensities = max_intensities[CNs > 0] ** 10
    print(max_intensities / sum(max_intensities))

    columns = np.zeros(image.shape)
    for loc in max_xys:
        x, y = loc * resolution
        columns[y][x] = 1
    
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, ax = plt.subplots(num=3)
    fig.colorbar(ax.pcolormesh(columns, cmap=cm.viridis, linewidths=0))
    ax.set_xlim((0, STEM_parameters['dimensions'][0] * 10))
    ax.set_ylim((0, STEM_parameters['dimensions'][1] * 10))    
        
    plt.show()
    print(len(min_intensities))
    import sys; sys.exit()
