import os
import math
import logging
import numpy as np
from scipy.signal import fftconvolve
from scipy.ndimage import sobel
from scipy.optimize import fmin

from ase.io import read

from structopt.tools import root, single_core, parallel
from structopt.tools.dictionaryobject import DictionaryObject
import gparameters

class STEM(object):
    """Calculates the chi^2 difference between a simulated and experimental image.
    In order to calculate a z-contrast image and chi^2 function the following
    parameters are required

    Parameters
    ----------
    HWHM : float
        The half-width half-maximum of the gaussian function used for the 
        point spread function.
    dimensions : list
        The x and y dimensions of STEM image.
    resolution : float
        The pixels per angstrom resolution
    """

    @single_core
    def __init__(self, parameters=None):
        if parameters is None:
            parameters = {}
        self.parameters = parameters
        self.parameters.setdefault('kwargs', {})
        self.parameters.kwargs.setdefault('zed', 1)
        self.psf = None
        self.target = None
        self.phantom = True

        # If running within StructOpt, create directory for saving files
        # and faster loading of PSF and target data
        path = os.path.join(gparameters.logging.path, 'fitness/STEM')
        path = os.path.join(path, 'rank-{}'.format(gparameters.mpi.rank))
        self.path = path
        os.makedirs(self.path, exist_ok=True)

    def calculate_fitness(self, individual):
        """Calculates the fitness of an individual with respect to a target
        image. Normalize this fitness by the number of atoms."""

        if self.psf is None:
            self.generate_psf()
        if self.target is None:
            self.generate_target()

        image = self.get_image(individual)
        image, x_shift, y_shift = self.cross_correlate(image)

        chi = image - self.target
        chi = self.normalize(chi, individual)

        return chi

    def cross_correlate(self, image):
        convolution = fftconvolve(self.target, image[::-1, ::-1], mode='full')
        y_max, x_max = np.unravel_index(np.argmax(convolution), convolution.shape)
        x_shift = x_max - image.shape[1] + 1
        y_shift = y_max - image.shape[0] + 1
        image = np.roll(image, x_shift, axis=1)
        image = np.roll(image, y_shift, axis=0)

        return image, x_shift, y_shift

    def normalize(self, chi, individual):
        if 'normalize' not in self.parameters:
            return chi

        norms = self.parameters['normalize']

        if 'SSE' in norms and norms['SSE']:
            chi = np.sum(np.square(chi)) ** 0.5
            if 'nprotons' in norms and norms['nprotons']:
                chi /= sum(individual.get_atomic_numbers()) ** 0.5
        else:
            chi = np.sum(np.absolute(chi))
            if 'nprotons' in norms and norms['nprotons']:
                chi /= sum(individual.get_atomic_numbers())

        return chi

    def get_Z_diff(self, individual):
        """Returns the difference between the individual and target's
        estimated atomic number count"""

        if self.target is None:
            self.generate_target()

        image = self.get_image(individual)
        zed = self.parameters.kwargs['zed']
        image_Z_tot = np.sum(image ** (1/zed))
        if self.phantom == True:
            target_Z_tot = np.sum(self.target ** (1/zed))
        else:
            return NotImplementedError("Real STEM image not implemented yet")
        
        Z_diff = image_Z_tot - target_Z_tot

        return Z_diff

    def get_linear_convolution(self, individual):
        """Calculate linear convoluted potential of an individual"""

        r = self.parameters.kwargs['resolution']
        zed = self.parameters.kwargs['zed']
        xmax, ymax = self.parameters.kwargs['dimensions']
        if isinstance(xmax, float):
            nx = int(xmax * r)
            ny = int(ymax * r)
            dx = xmax/nx
            dy = ymax/ny
        elif isinstance(xmax, int):
            nx = xmax
            ny = ymax
            xmax = nx / r
            ymax = ny / r
            dx = xmax / nx
            dy = ymax / ny

        ax, ay, az = individual.get_positions().T

        # Assign atom to the bottom left of the grid point
        ix, iy = np.floor(ax / dx), np.floor(ay / dy)
        
        # Apply periodic boundary conditions, considering all
        # corners of each pixel
        iax, ibx = np.fmod(ix, nx), np.fmod(ix + 1, nx)
        iay, iby = np.fmod(iy, ny), np.fmod(iy + 1, ny)

        # Array of fraction of atoms at the left and bottom of the pixel
        fax = 1 - np.fmod(ax / dx, 1)
        fay = 1 - np.fmod(ay / dy, 1)

        # Add potentials to grid. Split up the potential into
        # fractions on the pixel
        Zatom = individual.get_atomic_numbers()
        V1 = fax * fay * Zatom ** zed
        V2 = (1 - fax) * fay * Zatom ** zed
        V3 = fax * (1 - fay) * Zatom ** zed
        V4 = (1 - fax) * (1 - fay) * Zatom ** zed

        V = np.zeros([nx, ny])
        for j in range(len(individual)):
            V[int(iax[j]), int(iay[j])] += V1[j]
            V[int(ibx[j]), int(iay[j])] += V2[j]
            V[int(iax[j]), int(iby[j])] += V3[j]
            V[int(ibx[j]), int(iby[j])] += V4[j]

        return V

    def generate_psf(self):
        """Generates a psf array built from a gaussian function. The relevant 
        parameters specified in the parameters dictionary are below."""

        # We do not want to generate the psf if it has already been saved
        if (self.path is not None
            and os.path.isfile(os.path.join(self.path, 'psf.npy'))):
            with open(os.path.join(self.path, 'psf.npy'), "rb") as npy:
                self.psf = np.load(npy)
            return

        HWHM = self.parameters.kwargs['HWHM']
        r = self.parameters.kwargs['resolution']
        a, b = self.parameters.kwargs['dimensions']

        if isinstance(a, float):
            N_a = int(a * r)
            N_b = int(b * r)
        else:
            N_a = a
            N_b = b

        k_a_min = -0.5 * r
        k_a_max = (0.5 - 1.0/N_a) * r
        k_b_min = -0.5 * r
        k_b_max = (0.5 - 1.0/N_b) * r

        k_a = np.linspace(k_a_min, k_a_max, N_a)
        k_b = np.linspace(k_b_min, k_b_max, N_b)

        Mk_a, Mk_b = np.meshgrid(k_a, k_b)
        Mksq = Mk_a ** 2 + Mk_b ** 2

        d_k =  (2 * np.log(2)) ** 0.5 / (HWHM * 2 * np.pi)

        psf = np.exp(-Mksq / (2 * d_k ** 2))

        self.psf = psf

        # Try saving the PSF for future calculations
        if self.path is not None:
            np.save(os.path.join(self.path, 'psf'), self.psf)

        return

    def generate_target(self):
        """Generates the target STEM image from the parameters. The bulk of
        this code generates the target from an atoms object. """

        # Load the target from a file
        if (self.path is not None
            and os.path.isfile(os.path.join(self.path, 'target.npy'))):
            with open(os.path.join(self.path, 'target.npy'), "rb") as npy:
                self.target = np.load(npy)            
            return

        if self.psf is None:
            self.generate_psf()

        if not self.parameters.kwargs['target'].endswith('.xyz'):
            self.target = self.read_target(self.parameters.kwargs['target'])
            self.phantom = False
        else:
            atoms = read(self.parameters.kwargs['target'])
            self.target = self.get_image(atoms)
            self.phantom = True

        # Try saving the target for future calculations
        if self.path is not None:
            np.save(os.path.join(self.path, 'target'), self.target)

        return

    def read_target(self, path):
        if path.endswith('.npy'):
            target = np.load(path)

        return target

    def get_image(self, individual):
        """Calculates the z-contrasted STEM image of an individual"""

        if self.psf is None:
            self.generate_psf()

        psf = self.psf
        V = self.get_linear_convolution(individual)

        ft_psf = np.fft.fftshift(psf)
        ft_V = np.fft.fft2(V).T

        image = np.fft.ifft2(ft_psf * ft_V, axes=(0, 1)).real

        if 'multislice' in self.parameters.kwargs:
            image = self.get_multislice(image, self.parameters.kwargs['multislice'])

        return image

    def get_multislice(self, image, multislice_params):
        """Converts pixel by pixel""" 
        coeffs = multislice_params['coeffs']
        plot_type = multislice_params['plot_type']

        if 'fit_resolution' not in multislice_params:
            scale = 50 / 3.9231
        else:
            scale = multislice_params['fit_resolution']

        image *= (self.parameters.kwargs['resolution'] / scale) ** 2

        if plot_type == 'log':
            image = np.log(image + 1)
            image = np.nan_to_num(np.poly1d(coeffs)(image))
        elif plot_type == 'log-d':
            log_image = np.log(image + 1)

            log_dx = sobel(log_image, axis=0, mode='wrap')
            log_dy = sobel(log_image, axis=1, mode='wrap')
            log_dxy = np.hypot(log_dx, log_dy)

            shape = image.shape
            log_image = log_image.flatten()
            log_dxy = log_dxy.flatten()

            A0 = np.ones(log_image.shape)
            A1 = log_image
            A2 = log_image ** 2
            A3 = log_image ** 3
            A4 = log_dxy

            A = np.vstack([A4, A3, A2, A1, A0]).T

            X = coeffs

            image = np.dot(A, X)
            image = np.reshape(image, shape)

        return image

