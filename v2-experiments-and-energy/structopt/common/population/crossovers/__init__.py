import functools
import random
from itertools import accumulate
from bisect import bisect

from structopt.tools import root, single_core, parallel, allgather
import gparameters

from .rotate import rotate

rotate.tag = 'Ro'


class Crossovers(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters

        # self.crossovers is a dictionary containing {function: probability} pairs
        self.crossovers = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters}

        # self.kwargs is a dictionary containing {function: kwargs} pairs
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters}

        self.total_probability = sum(self.crossovers.values())
        assert self.total_probability <= 1.0
        self.crossovers[None] = 1.0 - self.total_probability
        self.selected_crossover = None


    @single_core
    def select_crossover(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.crossovers.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_crossover = choices[bisect(cumdist, x)]

    @parallel
    def crossover(self, pairs):
        """ """
        ncores = gparameters.mpi.ncores
        rank = gparameters.mpi.rank

        # Assign which pairs to mate on which cores
        pairs_per_core = {rank: [] for rank in range(ncores)}
        for i, pair in enumerate(pairs):
            pairs_per_core[i % ncores].append(pair)

        # Perform the designated crossovers by rank
        children = []
        for individual1, individual2 in pairs_per_core[rank]:
            self.select_crossover()  # Choose a new crossover to perform for every pair
            if self.selected_crossover is not None:
                kwargs = self.kwargs[self.selected_crossover]
                child1, child2 = self._crossover(individual1, individual2, self.selected_crossover, kwargs)
                children.append(child1)
                children.append(child2)
            else:
                children.append(None)
                children.append(None)

        children_per_core = {r: [] for r in range(ncores)}
        all_children = []
        count = 0
        for i in range(ncores):
            if i == rank:
                for child in children:
                    all_children.append(child)
                    children_per_core[i].append(count)
                    count += 1
            else:
                for _ in pairs_per_core[i]:
                    all_children.append(None)
                    children_per_core[i].append(count)
                    count += 1
                    all_children.append(None)
                    children_per_core[i].append(count)
                    count += 1

        if gparameters.mpi.ncores > 1:
            all_children = allgather(all_children, children_per_core)

        count_nones = all_children.count(None)
        all_children = [child for child in all_children if child is not None]

        assert len(all_children) == len(pairs)*2 - count_nones

        return all_children

    @single_core
    def _crossover(self, individual1, individual2, crossfunction, crosskwargs):
        if crossfunction is None:
            raise ValueError("Tried to perform a crossover but the selected crossover was `None`.")
        print("Performing crossover {} on individuals {} and {}".format(crossfunction.__name__, individual1, individual2))
        child1, child2 = crossfunction(individual1, individual2, **crosskwargs)
        if child1 is not None:
            child1._fitted = False
            child1._relaxed = False
        if child2 is not None:
            child2._fitted = False
            child2._relaxed = False
        self.post_processing((individual1, individual2), (child1, child2))
        return child1, child2


    @single_core
    def post_processing(self, parent_pair, child_pair):
        parent1, parent2 = parent_pair
        child1, child2 = child_pair
        for child in child_pair:
            if child is not None:
                child.crossover_tag = 'c{tag}({parent1}+{parent2})'.format(parent1=parent1.id, parent2=parent2.id, tag=self.selected_crossover.tag)


    @staticmethod
    @functools.wraps(rotate)
    def rotate(individual1, individual2, conserve_composition=True):
        return rotate(individual1, individual2, conserve_composition)

