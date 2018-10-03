from structopt.common.crossmodule import get_avg_radii

def get_particle_radius(atomlist, fill_factor=0.74):
    """Returns an estimated nanoparticle radius given a
    concentration of atoms and void fraction of the sphere.
    Given an average sphere, this is given by the formula

    R_sphere = (n_tot / f)**(1.0/3.0) * R_atom

    where n_tot is the total number of atoms and f is
    the fill factor of the particle.

    Parameters
    ----------
    atomlist : list
        An N x M list where N is the number of unique atoms and M are
        the attributes of the atom. The first and second index of each
        N list must be the chemical symbol and number, respectively.

    fill_factor : float between 0.0 and 1.0
        Factor that determines the packing of spheres in the particle.

    Output
    ------
    out : float
        Radius of nanoparticle.
    """

    n_tot = sum([atom[1] for atom in atomlist])
    R_atom = get_avg_radii(atomlist)
    R_sphere = (n_tot / fill_factor)**(1.0/3.0) * R_atom

    return R_sphere
