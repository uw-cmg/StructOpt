import logging
import functools
import random
from itertools import accumulate
from bisect import bisect
from collections import defaultdict

import structopt
from structopt.tools import root, single_core, parallel

from .swap_positions import swap_positions
from .swap_species import swap_species
from .move_atoms import move_atoms
from .rotate_atoms import rotate_atoms
from .rotate_cluster import rotate_cluster
from .rotate_all import rotate_all
from .permutation import permutation
from .rattle import rattle

swap_positions.tag = 'SwPo'
swap_species.tag = 'SwSp'
move_atoms.tag = 'MoAt'
rotate_atoms.tag = 'RoAt'
rotate_cluster.tag = 'RoCl'
rotate_all.tag = 'RoAl'
permutation.tag = 'Pe'
rattle.tag = 'Rat'

NOT_MUTATIONS = ['preserve_best', 'keep_original', 'keep_original_best']

class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters

        # self.mutations is a dictionary containing {function: probability} pairs
        self.mutations = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters
                          if name not in NOT_MUTATIONS}

        #self.kwargs is a dictionary containing {function: kwargs} pairs
        self.kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters
                       if name not in NOT_MUTATIONS}

        total_probability = sum(self.mutations.values())
        assert total_probability <= 1.0
        self.mutations[None] = 1.0 - total_probability
        self.selected_mutation = None


    @single_core
    def select_mutation(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.mutations.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_mutation = choices[bisect(cumdist, x)]


    @single_core
    def mutate(self, individual):
        if self.selected_mutation is None:
            return

        logger = logging.getLogger("default")
        logger.info("Performing mutation {} on individual {}".format(self.selected_mutation.__name__, individual.id or getattr(individual, "mutated_from", None)))
        print("Performing mutation {} on individual {}".format(self.selected_mutation.__name__, individual.id or getattr(individual, "mutated_from", None)))

        kwargs = self.kwargs[self.selected_mutation]
        result = self.selected_mutation(individual, **kwargs)

        # If the mutation "failed" and therefore did not modify the individual, do not update the below attributes
        if result is False:
            print("FAILED performing mutation {} on individual {}".format\
                (self.selected_mutation.__name__, individual.id \
                or getattr(individual, "mutated_from", None)))
            logger.info("FAILED performing mutation {} on individual {}".format\
                (self.selected_mutation.__name__, individual.id \
                or getattr(individual, "mutated_from", None)))
            return individual

        individual._relaxed = False
        individual._fitted = False
        self.post_processing(individual)
        return


    @single_core
    def post_processing(self, individual):
        individual.mutation_tag = 'm{tag}({id})'.format(tag=self.selected_mutation.tag, id=getattr(individual, "mutated_from", "?"))


    @staticmethod
    @functools.wraps(swap_positions)
    def swap_positions(individual):
        return swap_positions(individual)


    @staticmethod
    @functools.wraps(swap_species)
    def swap_species(individual):
        return swap_species(individual)


    @staticmethod
    @functools.wraps(move_atoms)
    def move_atoms(individual):
        return move_atoms(individual)


    @staticmethod
    @functools.wraps(rotate_atoms)
    def rotate_atoms(individual):
        return rotate_atoms(individual)


    @staticmethod
    @functools.wraps(rotate_cluster)
    def rotate_cluster(individual):
        return rotate_cluster(individual)


    @staticmethod
    @functools.wraps(rotate_all)
    def rotate_all(individual, vector=None, angle=None, center=None):
        return rotate_all(individual, vector, angle, center)    


    @staticmethod
    @functools.wraps(permutation)
    def permutation(individual):
        return permutation(individual)

    @staticmethod
    @functools.wraps(rattle)
    def rattle(individual, stdev=0.5, x_avg_bond=True):
        return rattle(individual, stdev, x_avg_bond)
