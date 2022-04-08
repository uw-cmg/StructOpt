import numpy as np
from structopt.common.crossmodule import get_avg_radii

def rattle(individual, stdev=0.5, x_avg_bond=True):
    """Randomly displace all atoms in a random direction with a magnitude
    drawn from a gaussian distribution.
    
    Parameters
    ----------
    individual : Individual
        An individual
    stdev : float
        The standard deviation of the gaussian distribution to rattle
        all the atoms. If x_avg_bond is set to True, given as the fraction
        of the average bond length of the material.
    x_avg_bond : bool
        If True, the gaussian distributions standard deviation is 
        stdev * avg_bond_length. Note, this only applies to fcc, hcp,
        or bcc materials.
    """

    if x_avg_bond:
        stdev *= get_avg_radii(individual) * 2

    pos = individual.get_positions()
    individual.set_positions(pos + np.random.normal(scale=stdev, size=pos.shape))

    return
