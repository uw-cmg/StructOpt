import random


def diversify_module(individual1, individual2, module='LAMMPS', min_diff=0.0001):
    """This module eliminates individuals that have similar fitness values.
    The goal is to diversify the population and eliminate duplicate 
    individuals. It is often useful to include this in the optimizer along
    with another predator.

    Parameters
    ----------
    individual1 : structopt.common.individual.Individual
        The first individual
    individual2 : structopt.common.individual.Individual
        The second individual
    module : str
        The name of the module for evaluating similar fitnesses.
    min_diff : float
        The minimum difference between individual fitness values for judging
        they are the same.
    """

    fit1 = getattr(individual1, module)
    fit2 = getattr(individual2, module)
    if abs(fit1 - fit2) < min_diff:
        return True
    else:
        return False

