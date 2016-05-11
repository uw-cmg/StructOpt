import functools
import random
from itertools import accumulate
from bisect import bisect

import structopt


class Crossovers(object):
    """ """
    def __init__(self):
        self.parameters = structopt.parameters.crossovers
        self.crossovers = {getattr(self, name): prob for name, prob in self.parameters.options.items()}
        total_probability = sum(self.crossovers.values())
        self.crossovers[None] = 1.0 - total_probability
        self.selected_crossover = None

    def select_crossover(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.crossovers.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_crossover = choices[bisect(cumdist, x)]

    def crossover(self, population):
        children = []
        for individual1, individual2 in zip(population[::2], population[1::2]):
            self.select_crossover()
            if self.selected_crossover is not None:
                child1, child2 = self.selected_crossover(individual1, individual2)
                if child1 is not None:
                    children.append(child1)
                if child2 is not None:
                    children.append(child2)
        return children
