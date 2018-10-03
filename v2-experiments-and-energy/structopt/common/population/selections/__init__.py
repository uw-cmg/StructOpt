import functools
import random
from itertools import accumulate
from bisect import bisect

from .random_selection import random_selection
from .rank import rank
from .roulette import roulette
from .tournament import tournament
from .best import best
from structopt.tools import root, single_core, parallel


class Selections(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.selections = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}
        total_probability = sum(self.selections.values())
        assert total_probability <= 1.0
        self.selections[None] = 1.0 - total_probability
        self.selected_selection = None


    @single_core
    def select_selection(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.selections.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_selection = choices[bisect(cumdist, x)]


    @single_core
    def select(self, population):
        fits = [individual.fitness for individual in population]
        if self.selected_selection is None:
            return []
        kwargs = self.kwargs[self.selected_selection]
        pairs = self.selected_selection(population=population, fits=fits, **kwargs)
        self.post_processing(pairs)
        return pairs


    @single_core
    def post_processing(self, pairs):
        pass


    @staticmethod
    @functools.wraps(random_selection)
    def random_selection(population, fits):
        return random_selection(population, fits)

    @staticmethod
    @functools.wraps(rank)
    def rank(population, fits, p_min=None, unique_pairs=False, unique_parents=False):
        return rank(population, fits, p_min, unique_pairs, unique_parents)

    @staticmethod
    @functools.wraps(roulette)
    def roulette(population, fits, unique_pairs=False, unique_parents=False):
        return roulette(population, fits, unique_pairs, unique_parents)

    @staticmethod
    @functools.wraps(tournament)
    def tournament(population, fits, tournament_size=5, unique_pairs=False, unique_parents=False, keep_best=False):
        return tournament(population, fits, tournament_size, unique_pairs, unique_parents, keep_best)

    @staticmethod
    @functools.wraps(best)
    def best(population, fits):
        return best(population, fits)
