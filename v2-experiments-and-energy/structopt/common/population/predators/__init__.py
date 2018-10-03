import functools
import logging
import random
from itertools import accumulate
from bisect import bisect
import numpy as np

from structopt.tools import root, single_core, parallel
from .best import best
from .roulette import roulette
from .tournament import tournament
from .rank import rank
from .fuss import fuss
import gparameters


class Predators(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.predators = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}
        total_probability = sum(self.predators.values())
        assert total_probability <= 1.0
        self.predators[None] = 1.0 - total_probability
        self.selected_predator = None

    @single_core
    def select_predator(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.predators.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_predator = choices[bisect(cumdist, x)]

    @single_core
    def kill(self, population, nkeep, keep_best=True):
        """Removes some individuals from the population.

        Parameters
        ----------
        nkeep : int
            The number of individuals to keep.
        keep_best : bool
            If set to True, the best individual is always included in the
            following generation.

        Returns
        -------
            The individuals that were removed from the population.
        """
        if self.selected_predator is None or len(population) <= nkeep:
            return []

        fits = {individual.id: individual.fitness for individual in population}
        if keep_best:
            best_id = min(fits, key=fits.get)
            fits.pop(best_id)
            nkeep -= 1

        kwargs = self.kwargs[self.selected_predator]
        to_keep = self.selected_predator(fits=fits, nkeep=nkeep, **kwargs)
        if keep_best:
            try:
                to_keep.append(best_id)
            except AttributeError:
                to_keep = np.append(to_keep, best_id)

        killed = [individual for individual in population if individual.id not in to_keep]
        new_population = [population[id] for id in to_keep]
        population.replace(new_population)

        self.post_processing(killed)
        return killed

    @single_core
    def post_processing(self, killed):
        logger = logging.getLogger("output")
        logger.info("Generation {}: Killed: {}".format(gparameters.generation, killed))

    @staticmethod
    @functools.wraps(best)
    def best(fits, nkeep):
        return best(fits, nkeep)

    @staticmethod
    @functools.wraps(roulette)
    def roulette(fits, nkeep, T=None):
        return roulette(fits, nkeep, T)

    @staticmethod
    @functools.wraps(tournament)
    def tournament(fits, nkeep, tournament_size=5):
        return tournament(fits, nkeep, tournament_size)

    @staticmethod
    @functools.wraps(rank)
    def rank(fits, nkeep, p_min=None):
        return rank(fits, nkeep, p_min)

    @staticmethod
    @functools.wraps(fuss)
    def fuss(fits, nkeep, nbest=1, fusslimit=10):
        return fuss(fits, nkeep, nbest=1, fusslimit=10)

