import logging
import numpy as np
import scipy.ndimage.filters as filters
from scipy.ndimage import center_of_mass
from scipy.optimize import fmin, brute

import structopt.common.individual.fitnesses
from structopt.tools import root, single_core, parallel
from structopt.tools import rotation_matrix
from structopt.common.crossmodule import get_avg_radii, NeighborList

import gparameters

class STEM(structopt.common.individual.fitnesses.STEM):
    """Rotates the individual to obtain a better match with STEM image.
    This is done by taking one bright column in the bulk region and reading
    its nearest neighbors as a projection into the xy plane. A bulk atom close
    to the individual COM and its nearest neighbors in the individual is selected
    and rotated to optimize the matching of its projection in the xy plane
    with that of the STEM image.

    Parameters
    ---------
    HWHM : float
        The half-width half-maximum of the gaussian function used for the 
        point spread function.
    dimensions : list
        The x and y dimensions of STEM image.
    resolution : float
        The pixels per angstrom resolution
    gridsize : int
        Given a individual nearest neighbor unit to be optimized with the
        STEM image, the gridsize determines how fine the search for
        the rotation be. Tests indicate a gridsize of 10 is suitable.
    """

    def __init__(self, parameters=None):
        if parameters is None:
            parameters = {'kwargs': {}}
        parameters['kwargs'].setdefault('rotation_grid', 10)
        parameters['kwargs'].setdefault('rotation_iterations', 2)
        parameters['kwargs'].setdefault('surface_moves', 10)
        parameters['kwargs'].setdefault('filter_size', 1)

        super().__init__(parameters)

    def relax(self, individual):
        if self.psf is None:
            self.generate_psf()
        if self.target is None:
            self.generate_target()

        # Align the atom to produce optimum matching with STEM
        self.align(individual)

        current_fitness = self.calculate_fitness(individual)
        rank = gparameters.mpi.rank
        print("Relaxing individual {} on rank {} with STEM".format(individual.id, rank))

        # Relax the atom by rotating it
        steps = self.parameters.kwargs['rotation_grid']
        for i in range(self.parameters.kwargs['rotation_iterations']):
            bonds = self.get_bulk_bonds(individual)
            projection = self.get_STEM_projection(individual)
            solution = brute(self.epsilon,
                             args=(bonds, projection),
                             ranges=(slice(0, np.pi, np.pi/steps),
                                     slice(-1, 1, 2/steps),
                                     slice(0, 2*np.pi, 2*np.pi/steps)),
                             finish=fmin)

            phi, costheta, a = solution
            theta = np.arccos(costheta)

            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)

            individual.rotate([x, y, z], a, center='COP')
            new_fitness = self.calculate_fitness(individual)
            if new_fitness > current_fitness:
                individual.rotate([x, y, z], -a, center='COP')
            else:
                self.iterations = i
                break

        # Align the atom to produce optimum matching with STEM
        self.align(individual)

        if hasattr(individual, 'id'):
            print("Finished relaxing individual {} on rank {} with STEM".format(individual.id, rank))

        return

    @staticmethod
    def epsilon(rotation, bonds, projection):
        """Calculates the difference in the projected xy coordinates
        of a rotated bonding center of that of the individual
        and of the STEM image

        Parameters
        ----------
        rotation : list [phi, cos(theta), a]
            Contains phi (rotation around z-axis) and costheta (rotation 
            along the z-axis) and the angle to rotate the atoms object
        bonds : list of [vx, vy, vz] vectors
            The bonds around an atomic center in the bulk of the individual.
        projection : list of [vx, vy] vectors
            The bonds around an atomic center from the STEM image

        Output
        ------
        out : float
            The sum squared difference between projected [vx, vy] vectors
            after performing rotation.
        """

        phi, costheta, a = rotation
        theta = np.arccos(costheta)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        rotate = rotation_matrix(np.array([x, y, z]), a)
        bonds = np.asarray([[np.dot(rotate, v)[:2] for v in bonds]])
        projection = np.asarray([projection])
        dists = np.linalg.norm(np.transpose(bonds, (1, 0, 2)) - projection, axis=2)
        total_error = np.sum(np.square(np.amin(dists, axis=1)))

        return total_error

    def get_bulk_bonds(self, individual):
        """Chooses a random atom to orient around. The probability
        of choosing the atom is proportional to the atom's distance
        from the center of mass"""

        # Get a bulk atom near the center of the particle
        pos = individual.get_positions()
        com = individual.get_center_of_mass()
        dists_from_com = np.linalg.norm(pos - com, axis=1)
        prob = dists_from_com / sum(dists_from_com)
        bulk_atom_index = np.random.choice(list(range(len(individual))), p=prob)
        bulk_atom_pos = pos[bulk_atom_index]

        # Get the neighbors of the bulk atom
        NNs = NeighborList(individual)
        bonds = np.asarray([pos[i] for i in NNs[bulk_atom_index]]) - bulk_atom_pos

        return bonds

    def get_STEM_projection(self, individual, test=False):
        """Gets the projection of bonds from the brighest STEM image"""

        if self.target is None:
            self.generate_target()

        parameters = self.parameters
        target = self.target

        # Get a cutoff between maximum points in the STEM image based
        # on nearest neighbor distances
        cutoff = get_avg_radii(individual) * 2 * 1.1
        size = cutoff * parameters.kwargs['resolution'] * parameters.kwargs['filter_size']

        # Get a list of xy positions from analyzing local maxima in STEM image
        # as well as the position of a spot near the center of mass
        data_max = filters.maximum_filter(target, size=size)
        maxima = ((target == data_max) & (target > 0.1)) # Filter out low maxima
        com = np.asarray(center_of_mass(target)[::-1]) / parameters.kwargs['resolution']
        pos = np.argwhere(maxima)[:,::-1] / parameters.kwargs['resolution']
        dists_from_com = np.linalg.norm(pos - com, axis=1)
        prob = dists_from_com / sum(dists_from_com)
        bulk_atom_index = np.random.choice(list(range(len(dists_from_com))), p=prob)
        bulk_atom_pos = pos[bulk_atom_index]
        vecs = pos - bulk_atom_pos
        dists = np.linalg.norm(vecs, axis=1)
        vecs = vecs[((dists < cutoff) & (dists > 0))]

        if test:
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm

            fig, ax = plt.subplots(num=1)
            fig.colorbar(ax.pcolormesh(target, cmap=cm.viridis, linewidths=0))
            ax.set_xlim((0, parameters.kwargs['dimensions'][0] * parameters.kwargs['resolution']))
            ax.set_ylim((0, parameters.kwargs['dimensions'][1] * parameters.kwargs['resolution']))

            fig, ax = plt.subplots(num=2)
            fig.colorbar(ax.pcolormesh(data_max, cmap=cm.viridis, linewidths=0))
            ax.set_xlim((0, parameters.kwargs['dimensions'][0] * parameters.kwargs['resolution']))
            ax.set_ylim((0, parameters.kwargs['dimensions'][1] * parameters.kwargs['resolution']))

            fig, ax = plt.subplots(num=3)
            fig.colorbar(ax.pcolormesh(maxima, cmap=cm.viridis, linewidths=0))
            ax.set_xlim((0, parameters.kwargs['dimensions'][0] * parameters.kwargs['resolution']))
            ax.set_ylim((0, parameters.kwargs['dimensions'][1] * parameters.kwargs['resolution']))

            plt.show()

        return vecs

    def align(self, atoms):
        x, y = fmin(self.chi2, [0, 0], args=(atoms, self), disp=False)
        atoms.translate([x, y, 0])
        return

    @staticmethod
    def chi2(shift, atoms, module):
        x, y = shift
        atoms = atoms.copy()
        atoms.translate([x, y, 0])
        image = module.get_image(atoms)
        return np.sum(np.square(module.target - image)) ** 0.5
