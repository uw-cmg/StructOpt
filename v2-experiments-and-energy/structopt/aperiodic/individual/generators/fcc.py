import random
from math import cos, sin
import numpy as np
np.seterr(all='ignore')

from ase import Atoms
from structopt.tools import random_three_vector

def fcc(atomlist, cell, a, shape=[1, 1, 1], orientation=None, size=21,
        roundness=0.5, alpha=10, v=None, angle=None):
    """Generates a fcc nanoparticle of varying shape and orientation.
    For multi-component particles, elements are randomly distributed.
    
    Parameters
    ----------
    atomlist : list
        A list of [sym, n] pairs where sym is the chemical symbol
        and n is the number of of sym's to include in the individual
    cell : list
        A 3 element list that defines the x, y, and z dimensions of the
        simulation box
    a : float
        fcc lattice constant
    shape : list
        The ratio of the x, y, and z dimensions of the particle.
    orientation : str
        The facet that is parallel to the xy plane. This is useful for
        LAMMPS+STEM calculations where one already knows the orientation.
        If None, a random orientation is chosen.
    size : int
        Size of the fcc grid for building the nanoparticle. For example,
        a size of 21 means 21 x 21 x 21 supercell of fcc primitive cells
        will be used. Note, if the size is too small, the function will
        automatically expand the cell to fit the particle.
    roundness : float
        Determines the "roundness" of the particle. A more round particle
        will have a smaller surface area to volume ratio, but more
        undercoordinated surface sites. A less round particle will take 
        more and more the form of an octahedron.
    alpha : float
        Parameter for determining how defective the particle will be. 
        Higher alpha means less defective particle.
    v : list
        Used for a custom orientation. V is the vector in which to rotate
        the particle. Requires angle parameters to be entered. All rotations
        are done with respect to the 100 plane.
    angle : float
        Angle, in radians, to rotate atoms around vector v. All rotations are 
        done with respect to the 100 plane.
    """

    grid = np.zeros((size, size, size), dtype=np.int)

    # Start the atom
    middle = size // 2
    center = [middle, middle, middle]
    grid[middle, middle, middle] = 1

    n = sum([atom[1] for atom in atomlist])

    v, angle = get_vector_angle(orientation, v, angle)

    dists_array = get_norm_dists(grid, center, shape, a, v=v, angle=angle)

    for i in range(n - 1):

        # Get an array of 0 - 1 fitnesses that depend on the vacancy coord
        # Fitness close to 1 corresponds to high fitness (high CN)
        coords = get_coordination_numbers(grid)
        vac_coords = (1 - grid) * coords

        # Do not consider coordination numbers less than 3
        if np.max(vac_coords) >= 3:
            np.place(vac_coords, vac_coords < 3, 0)
        max_coord = np.max(vac_coords)
        min_coord = np.min(vac_coords[np.nonzero(vac_coords)])

        if max_coord != min_coord:
            max_array = max_coord * vac_coords.astype(bool).astype(float)
            min_array = min_coord * vac_coords.astype(bool).astype(float)
            coord_fit = np.nan_to_num((vac_coords-min_array) / (max_array-min_array))
        else:
            coord_fit = vac_coords // np.max(vac_coords)        

        # Get an array 0 - 1 fitness that depend on the distance
        # From the center of the atom. Fitness close to 1 corresponds
        # to high fitness (close to center)
        vac_grid = (1 - grid) * coords.astype(bool).astype(int)
        vac_dists = dists_array * vac_coords

        max_dist = np.max(vac_dists)
        min_dist = np.min(vac_dists[np.nonzero(vac_dists)])
        if max_dist != min_dist:
            dists_fit = min_dist / vac_dists
            dists_fit[dists_fit == np.inf] = 0
        else:
            dists_fit = vac_dists.astype(bool).astype(float)

        # Combine fitnesses using a roundness parameter. Weight highest
        # fitnesses according to a polynomial distribution
        total_fit = roundness*dists_fit + (1 - roundness)*coord_fit
        max_fit = np.nanmax(total_fit)
        min_fit = np.nanmin(total_fit)
        if max_fit != min_fit:
            total_fit = ((total_fit - min_fit) / (max_fit - min_fit))**(alpha)
        else:
            total_fit = total_fit.astype(bool).astype(float)

        # Normalize the probabilities and flatten the array
        sum_fit = np.nansum(total_fit)
        grow_prob = total_fit.flatten()
        grow_prob = (grow_prob / sum_fit)

        # Choose vacancy with the fitness and back calculate the [z, y, x]
        # index using the flattened index. Finally add the atom
        add_ind = np.random.choice(np.arange(len(grow_prob))[grow_prob > 0],
                                   p=grow_prob[grow_prob > 0])
        nz, ny, nx = np.shape(grid)
        add_ind_x = add_ind % (nx * ny) % nx
        add_ind_y = (add_ind % (nx * ny) - add_ind_x) / nx
        add_ind_z = (add_ind - add_ind_x - add_ind_y * nx) / (nx * ny)
        add_ind = np.asarray((add_ind_z, add_ind_y, add_ind_x), dtype=int)
        grid[add_ind[0], add_ind[1], add_ind[2]] = 1

        # Expand the grid if it's not large enough
        if (0 in add_ind or size - 1 in add_ind):
            for axis in range(3):
                expand = list(np.shape(grid))
                expand[axis] = 2
                expand = np.zeros(expand, dtype=int)
                grid = np.concatenate((grid, expand), axis=axis)
                grid = np.roll(grid, 1, axis=axis)
            
            center = [i + 1 for i in center]
            dists_array = get_norm_dists(grid, center, shape, a, v=v, angle=angle)

    return get_atoms(grid, atomlist, cell, a, v=v, angle=angle)

def get_coordination_numbers(grid):
    '''Returns the coordination number of every position in a fcc grid'''

    shape = np.add([2, 2, 2], np.shape(grid))
    neighbors = np.zeros(shape, dtype=np.int)
    for sides in [((1, 1), (1, 1), (2, 0)),
                  ((1, 1), (1, 1), (0, 2)),
                  ((1, 1), (2, 0), (1, 1)),
                  ((1, 1), (0, 2), (1, 1)),
                  ((2, 0), (1, 1), (1, 1)),
                  ((0, 2), (1, 1), (1, 1)),
                  ((1, 1), (2, 0), (0, 2)),
                  ((2, 0), (1, 1), (0, 2)),
                  ((1, 1), (0, 2), (2, 0)),
                  ((2, 0), (0, 2), (1, 1)),
                  ((0, 2), (2, 0), (1, 1)),
                  ((0, 2), (1, 1), (2, 0))]:
        neighbors = np.add(neighbors, np.pad(grid, sides, mode='constant'))

    return neighbors[1:-1, 1:-1, 1:-1]

def get_atoms(grid, atomlist, cell, a, v=[1, 0, 0], angle=0.0):
    '''Returns an atoms object from a grid'''

    a1 = a/2.0 * np.array([1, 1, 0])
    a2 = a/2.0 * np.array([1, 0, 1])
    a3 = a/2.0 * np.array([0, 1, 1])
    temp_cell = np.array([a1, a2, a3])

    if np.shape(v) == (3,) and not hasattr(angle, '__iter__'):
        vs = np.expand_dims(v, 0)
        angles = [angle]
    else:
        vs = v
        angles = angle

    for v, angle in zip(vs, angles):
        v = np.asarray(v, dtype=float)
        v /= np.linalg.norm(v)

        c = cos(angle)
        s = sin(angle)

        temp_cell[:] = (c * temp_cell -
                        np.cross(temp_cell, s * v) +
                        np.outer(np.dot(temp_cell, v), (1.0 - c) * v))
    
    # Loop through each lattice point and add the atom
    scaled_positions = []
    size = np.shape(grid)
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                if grid[i][j][k] == 1:
                    scaled_positions.append([i, j, k])

    # Add the atoms to the atoms object
    size = len(scaled_positions)
    chemical_symbols = []
    for i in atomlist:
        chemical_symbols += [i[0]] * i[1]
    random.shuffle(chemical_symbols)        
    atoms = Atoms(chemical_symbols)
    atoms.set_cell(temp_cell)
    atoms.set_scaled_positions(scaled_positions)
    atoms.set_cell(cell)
    atoms.center()

    # Make the cell the size of the grid
    # A = np.array([np.shape(grid)])
    # cell = np.multiply(A.T, np.array([a1, a2, a3]))
    # atoms.set_cell(cell)

    return atoms

def get_norm_dists(grid, center, shape, a, v=[0, 0, 1], angle=0.0):
    '''Returns a grid of normalized distances from the center'''

    a1 = a/2.0 * np.array([1, 1, 0])
    a2 = a/2.0 * np.array([1, 0, 1])
    a3 = a/2.0 * np.array([0, 1, 1])
    cell = np.array([a1, a2, a3])

    if np.shape(v) == (3,) and not hasattr(angle, '__iter__'):
        vs = np.expand_dims(v, 0)
        angles = [angle]
    else:
        vs = v
        angles = angle

    for v, angle in zip(vs, angles):
        v = np.asarray(v, dtype=float)
        v /= np.linalg.norm(v)

        c = cos(angle)
        s = sin(angle)

        cell[:] = (c * cell -
                   np.cross(cell, s * v) +
                   np.outer(np.dot(cell, v), (1.0 - c) * v))

    # Get the coordinates of each grid point in real space
    size = np.shape(grid)
    inds = np.array([i for i in np.ndindex(size[0], size[1], size[2])])
    coords = np.dot(inds, cell) - np.dot(center, cell)

    # Now reweight these x,y,z distances by the dimensions of the box
    # These will be normalized to the x coordinate
    x_over_y = float(shape[0]) / float(shape[1])
    x_over_z = float(shape[0]) / float(shape[2])
    coords *= np.array([[1, x_over_y, x_over_z]])

    dists = np.linalg.norm(coords, axis=1)
    dists = np.reshape(dists, size)

    return dists

def get_vector_angle(orientation=None, v=None, angle=None):
    if (np.shape(v) == (3,) and type(a) in [float, int]):
        v = np.asarray(v, dtype=float)
        v /= np.linalg.norm(v)
        return v, angle

    elif orientation is None:
        angle = np.random.uniform(0,np.pi*2)
        v = random_three_vector()
        return v, angle

    elif orientation == '100':
        orientation_angle = 0.0
        orientation_v = np.array([1, 0, 0])
    elif orientation == '110':
        orientation_angle = np.pi / 4
        orientation_v = np.array([0, 1, 0])
    elif orientation == '111':
        orientation_angle = np.arcsin(1.0 / 3.0 ** 0.5)
        orientation_v = np.array([-1.0 / 2.0 ** 0.5, 1.0 / 2.0 ** 0.5, 0])
    else:
        raise NotImplementedError('Orientation not implemented')

    # Sometimes we want another rotation after an orientation
    if v is None and angle is not None:
        v = [orientation_v, np.array([0, 0, 1])]
        angle = [orientation_angle, np.asarray(angle)]
    else:
        v = orientation_v
        angle = orientation_angle        

    return v, angle
