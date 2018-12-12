import random

import numpy as np
from scipy.ndimage import filters
from ase import Atom, Atoms

from structopt.common.crossmodule import CoordinationNumbers
from structopt.common.crossmodule import get_avg_radii
from structopt.common.individual.fitnesses import STEM

def add_atom_STEM(individual, STEM_parameters, add_prob=None, permute=0.5, 
                  filter_size=1, column_cutoff=0.2, surf_cutoff=0.5, min_cutoff=0.5):

    """Moves surface atoms around based on the difference in the target
    and individual STEM image

    Parameters
    ----------
    STEM_parameters : dict
        Parameters for the STEM calculation. Ideally should be the same as the ones
        used for the STEM fitness/relaxation
    move_CN : int
        The maximum coordination number considered as a surface atom to move.    
    surf_CN : int
        The maximum coordination number considered as a surface atom to move an
        atom to
    filter_size : float
        Filter size for choosing local maximum in the picture. Filter size is equal
        to average_bond_length * resolution * filter_size.
    surf_cutoff : float
        The search radius for selecting a surface site near a low intensity point
        Defaults to the average bond distance
    """

    module = STEM({'kwargs': STEM_parameters})
    module.generate_target()
    target = module.target

    image, x_shift, y_shift = module.cross_correlate(module.get_image(individual))
    contrast = image - target
    min_min = np.min(contrast)

    # Determine filter size for locating local minimum
    avg_bond_length = get_avg_radii(individual) * 2
    cutoff = avg_bond_length * 1.1
    surf_cutoff *= cutoff
    column_cutoff *= cutoff
    resolution = module.parameters['kwargs']['resolution']
    size = cutoff * resolution * filter_size

    ###################################
    ## Code for testing the contrast ##
    ###################################
    # import matplotlib.pyplot as plt
    # import matplotlib.cm as cm
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(contrast, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0]*STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1]*STEM_parameters['resolution']))
    # plt.show()
    # import sys; sys.exit()

    # Get xy coordinates of the minimum intensities
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
    # fig, ax = plt.subplots()
    # fig.colorbar(ax.pcolormesh(minima, cmap=cm.viridis, linewidths=0))
    # ax.set_xlim((0, STEM_parameters['dimensions'][0]*STEM_parameters['resolution']))
    # ax.set_ylim((0, STEM_parameters['dimensions'][1]*STEM_parameters['resolution']))
    # plt.show()
    # print(len(min_intensities))
    # import sys; sys.exit()

    # Randomly choose local maxima and minima locations from contrast weighted
    # by their intensity
    low_xy_index = np.random.choice(np.arange(len(min_xys)), p=min_intensities)
    low_xy = min_xys[low_xy_index]

    # Get indices of atoms considered to be moved and sites to move to
    # Organize atoms into columns
    pos = individual.get_positions()
    xys = np.expand_dims(pos[:, :2], 0)
    dists = np.linalg.norm(xys - np.transpose(xys, (1, 0, 2)), axis=2)

    NNs = np.sort(np.argwhere(dists < column_cutoff))
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

        # Get the average bond length of each column for adding
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

    surf_xys = surf_positions[:,:2]

    dists_surf_xys = np.linalg.norm(surf_xys - low_xy, axis=1)
    indices_surf_xys = [i for i, d in enumerate(dists_surf_xys) if d < surf_cutoff]

    ########################
    ## Test atoms to move ##
    ########################
    # from ase import Atom
    # from ase.visualize import view
    # for i in indices_surf_xys:
    #     individual.append(Atom('Mo', surf_positions[i]))
    # view(individual)
    # import sys; sys.exit()

    # Pick the surface site to add the atom to
    if len(indices_surf_xys) == 0:
        surf_index = np.argmin(dists_surf_xys)
    else:
        surf_index = random.choice(indices_surf_xys)
    new_position = surf_positions[surf_index]
    new_xy = new_position[:2]

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

    # Now maybe switch the surface element with a bulk element
    if random.random() < permute:
        individual.append(Atom(element, new_position))
        return

    xys = individual.get_positions()[:,:2]
    syms = individual.get_chemical_symbols()
    dists_xys = np.linalg.norm(xys - new_xy, axis=1)
    indices_to_switch = [i for i, d in enumerate(dists_xys) if d < column_cutoff and syms[i] != element]

    if len(indices_to_switch) == 0:
        individual.append(Atom(element, new_position))
        return

    index_to_switch = random.choice(indices_to_switch)
    element_to_switch = syms[index_to_switch]
    individual[index_to_switch].symbol = element
    individual.extend(Atoms([Atom(element_to_switch, new_position)]))

    return
