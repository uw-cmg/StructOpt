import numpy as np
from scipy.signal import fftconvolve
from ase.cluster.octahedron import Octahedron
from ase.visualize import view

from structopt.cluster.individual.generators import fcc
from structopt.common.crossmodule.analysis import get_avg_radii

def get_chi2(atoms1, atoms2, cutoff=0.8, r=2.0, HWHM=0.4):
    """Calculates the chi2, which is the difference in positions
    between atoms1 and atoms2"""

    # Calculate the offset to apply to atoms1 to maximize matching
    # to atoms2
    offset = get_offset(atoms1, atoms2, r=r, HWHM=HWHM)
    atoms1.translate(offset)

    # Calculate the location difference of each atom1 atom with atom2 atom
    cutoff *= get_avg_radii(atoms1) * 2
    pos1 = np.expand_dims(atoms1.get_positions(), 0)
    pos2 = np.expand_dims(atoms2.get_positions(), 0)
    dists = np.linalg.norm(pos1 - np.transpose(pos2, (1, 0, 2)), axis=2)

    nn_dists = np.min(dists, axis=1)
    x_fp = np.sum(nn_dists > cutoff)
    x_fn = np.sum(np.min(dists, axis=0) > cutoff)
    chi2 = nn_dists[nn_dists <= cutoff]

    return x_fp, x_fn, chi2

def get_chi2_column(atoms1, atoms2, cutoff=0.2, r=2.0, HWHM=0.4):
    """Calculates the error per column of atoms in the z-direction"""

    cutoff *= get_avg_radii(atoms1) * 2

    # Calculate the offset to apply to atoms1 to maximize matching
    # to atoms2
    offset = get_offset(atoms1, atoms2, r=r, HWHM=HWHM)
    atoms1.translate(offset)

    # Group each atom in both atoms1 and atoms2 into columns
    xys1 = np.expand_dims(atoms1.get_positions()[:, :2], 0)
    xys2 = np.expand_dims(atoms2.get_positions()[:, :2], 0)
    dists1 = np.linalg.norm(xys1 - np.transpose(xys1, (1, 0, 2)), axis=2)
    dists2 = np.linalg.norm(xys2 - np.transpose(xys2, (1, 0, 2)), axis=2)

    NNs1 = np.sort(np.argwhere(dists1 < cutoff))
    column_indices1 = []
    atoms_to_be_sorted1 = list(range(len(atoms1)))
    while len(atoms_to_be_sorted1) > 0:
        i = atoms_to_be_sorted1[0]
        same_column_indices = np.unique(NNs1[NNs1[:,0] == i])
        column_indices1.append(same_column_indices)
        for j in reversed(sorted(same_column_indices)):
            i_del = atoms_to_be_sorted1.index(j)
            atoms_to_be_sorted1.pop(i_del)
            NNs1 = NNs1[NNs1[:,0] != j]
            NNs1 = NNs1[NNs1[:,1] != j]

    NNs2 = np.sort(np.argwhere(dists2 < cutoff))
    column_indices2 = []
    atoms_to_be_sorted2 = list(range(len(atoms2)))
    while len(atoms_to_be_sorted2) > 0:
        i = atoms_to_be_sorted2[0]
        same_column_indices = np.unique(NNs2[NNs2[:,0] == i])
        column_indices2.append(same_column_indices)
        for j in reversed(sorted(same_column_indices)):
            i_del = atoms_to_be_sorted2.index(j)
            atoms_to_be_sorted2.pop(i_del)
            NNs2 = NNs2[NNs2[:,0] != j]
            NNs2 = NNs2[NNs2[:,1] != j]            

    # Find the average xy coordinates of each column
    column_xys1 = np.asarray([np.average(np.squeeze(xys1, axis=0)[i], axis=0) for i in column_indices1])
    column_xys2 = np.asarray([np.average(np.squeeze(xys2, axis=0)[i], axis=0) for i in column_indices2])

    # Find matching column locations in atoms1 and atoms2
    dists = np.linalg.norm(np.expand_dims(column_xys1, 0) - np.transpose(np.expand_dims(column_xys2, 0), (1, 0, 2)), axis=2)
    paired_atoms2 = [np.argmin(dist) if np.min(dist) < cutoff else False for dist in dists]
    paired_atoms1 = [np.argmin(dist) if np.min(dist) < cutoff else False for dist in dists.T]

    n_fn = np.count_nonzero(np.array(paired_atoms2) == False)
    n_fp = np.count_nonzero(np.array(paired_atoms1) == False)
    column_pairs = sorted([[j, i] for i, j in enumerate(paired_atoms2) if j is not False], key=lambda i: i[0])

    syms1 = np.asarray(atoms1.get_chemical_symbols())
    syms2 = np.asarray(atoms2.get_chemical_symbols())

    unique_syms = np.unique(syms2)
    chi2 = {sym: [] for sym in unique_syms}
    chi2['n'] = []
    
    for pair in column_pairs:
        col_indices1 = column_indices1[pair[0]]
        col_indices2 = column_indices2[pair[1]]
        col_syms1 = list(syms1[col_indices1])
        col_syms2 = list(syms2[col_indices2])
        chi2['n'].append(len(col_syms1) - len(col_syms2))
        for sym in unique_syms:
            chi2[sym].append(col_syms1.count(sym) - col_syms2.count(sym))

    return n_fn, n_fp, chi2

def get_offset(atoms1, atoms2, r=5.0, HWHM=0.4):
    """Gets the offset to apply to atoms1 to have its positions match atoms2"""

    cell1 = atoms1.get_cell()
    cell2 = atoms2.get_cell()

    # Make sure the cell is a box and they are the same
    assert (cell1.diagonal() * np.eye(3) == cell1).all()
    assert (cell2.diagonal() * np.eye(3) == cell2).all()
    assert (cell1 == cell2).all()

    # Load the 3d point spread function (psf)
    psf = get_3d_psf(cell1.diagonal(), r, HWHM)
    ft_psf = np.fft.fftshift(psf)

    # Get the gridded locations
    V1 = get_gridded_locations(cell1.diagonal(), r, atoms1)
    ft_V1 = np.fft.fftn(V1).T
    image1 = np.fft.ifftn(ft_psf * ft_V1)

    V2 = get_gridded_locations(cell2.diagonal(), r, atoms2)
    ft_V2 = np.fft.fftn(V2).T
    image2 = np.fft.ifftn(ft_psf * ft_V2)

    # Use cross-correlation to calculate the ideal offset
    correlation = fftconvolve(image2, image1[::-1, ::-1, ::-1], mode='full')
    z_max, y_max, x_max = np.unravel_index(np.argmax(correlation), correlation.shape)
    x_shift = (x_max - image1.shape[2] + 1) / r
    y_shift = (y_max - image1.shape[1] + 1) / r
    z_shift = (z_max - image1.shape[0] + 1) / r

    return (x_shift, y_shift, z_shift)

def get_3d_psf(dimensions, r, HWHM):
    """Generates a psf array built from a gaussian function. The relevant 
    parameters specified in the parameters dictionary are below."""

    a, b, c = dimensions

    N_a = int(a * r)
    N_b = int(b * r)
    N_c = int(c * r)

    k_a_min = -0.5 * r
    k_a_max = (0.5 - 1.0/N_a) * r
    k_b_min = -0.5 * r
    k_b_max = (0.5 - 1.0/N_b) * r
    k_c_min = -0.5 * r
    k_c_max = (0.5 - 1.0/N_c) * r

    k_a = np.linspace(k_a_min, k_a_max, N_a)
    k_b = np.linspace(k_b_min, k_b_max, N_b)
    k_c = np.linspace(k_c_min, k_c_max, N_c)

    Mk_a, Mk_b, Mk_c = np.meshgrid(k_a, k_b, k_c)
    Mksq = Mk_a ** 2 + Mk_b ** 2 + Mk_c ** 2

    d_k =  (2 * np.log(2)) ** 0.5 / (HWHM * 2 * np.pi)

    psf = np.exp(-Mksq / (2 * d_k ** 2))

    return psf

def get_gridded_locations(dimensions, r, individual):
    """Calculate linear convoluted potential of an individual"""

    xmax, ymax, zmax = dimensions
    nx = int(xmax * r)
    ny = int(ymax * r)
    nz = int(zmax * r)
    dx = xmax/nx
    dy = ymax/ny
    dz = zmax/nz

    ax, ay, az = individual.get_positions().T

    # Assign atom to the bottom left of the grid point
    ix, iy, iz = np.floor(ax / dx), np.floor(ay / dy), np.floor(az / dz)

    # Apply periodic boundary conditions, considering all
    # corners of each pixel
    iax, ibx = np.fmod(ix, nx), np.fmod(ix + 1, nx)
    iay, iby = np.fmod(iy, ny), np.fmod(iy + 1, ny)
    iaz, ibz = np.fmod(iz, nz), np.fmod(iz + 1, nz)

    # Array of fraction of atoms at the left and bottom of the pixel
    fax = 1 - np.fmod(ax / dx, 1)
    fay = 1 - np.fmod(ay / dy, 1)
    faz = 1 - np.fmod(az / dz, 1)

    # Add potentials to grid. Split up the potential into
    # fractions on the pixel
    Zatom = np.ones(len(individual))
    V1 = fax * fay * faz * Zatom
    V2 = (1 - fax) * fay * faz * Zatom
    V3 = fax * (1 - fay) * faz * Zatom
    V4 = (1 - fax) * (1 - fay) * faz * Zatom
    V5 = fax * fay * (1 - faz) * Zatom
    V6 = (1 - fax) * fay * (1 - faz) * Zatom
    V7 = fax * (1 - fay) * (1 - faz) * Zatom
    V8 = (1 - fax) * (1 - fay) * (1 - faz) * Zatom

    V = np.zeros([nx, ny, nz])
    for j in range(len(individual)):
        V[int(iax[j]), int(iay[j]), int(iaz[j])] += V1[j]
        V[int(ibx[j]), int(iay[j]), int(iaz[j])] += V2[j]
        V[int(iax[j]), int(iby[j]), int(iaz[j])] += V3[j]
        V[int(ibx[j]), int(iby[j]), int(iaz[j])] += V4[j]
        V[int(iax[j]), int(iay[j]), int(iaz[j])] += V5[j]
        V[int(ibx[j]), int(iay[j]), int(iaz[j])] += V6[j]
        V[int(iax[j]), int(iby[j]), int(iaz[j])] += V7[j]
        V[int(ibx[j]), int(iby[j]), int(iaz[j])] += V8[j]

    return V
