import functools
import structopt.common.population.crossovers

from .rotate import rotate

rotate.tag = 'Ro'

class Crossovers(structopt.common.population.crossovers.Crossovers):

    @staticmethod
    @functools.wraps(rotate)
    def rotate(individual1, individual2, center_at_atom=True, repair_composition=True):
        return rotate(individual1, individual2, center_at_atom, repair_composition)
