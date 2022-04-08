import functools
import random
from itertools import accumulate, combinations
from bisect import bisect
from mpi4py import MPI
from structopt.tools import root, single_core, parallel, disjoint_set_merge
from structopt.tools.parallel import allgather
import gparameters
from .all_close_atom_positions import all_close_atom_positions
from .diversify_module import diversify_module


class Fingerprinters(object):
    """ """
    kwargs = ['keep_best']

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.fingerprinters = {getattr(self, name): self.parameters[name]['probability'] for name in self.parameters if name not in self.kwargs}
        self.function_kwargs = {getattr(self, name): self.parameters[name]['kwargs'] for name in self.parameters if name not in self.kwargs}
        total_probability = sum(self.fingerprinters.values())
        assert total_probability <= 1.0
        self.fingerprinters[None] = 1.0 - total_probability
        self.selected_fingerprinter = None

    @single_core
    def select_fingerprinter(self):
        # Implementation from https://docs.python.org/3/library/random.html -- Ctrl+F "weights"
        choices, weights = zip(*self.fingerprinters.items())
        cumdist = list(accumulate(weights))
        x = random.random() * cumdist[-1]
        self.selected_fingerprinter = choices[bisect(cumdist, x)]

    @single_core
    def remove_duplicates(self, population, nkeep, keep_best=True):
        assert nkeep <= len(population)
        if nkeep == len(population):
            return

        # Apply the selected fingerprinter to all pairs of individuals
        # The below function parallelizes the calculations, but returns the gathered results
        # So at this point, all cores are working with the same data and will run the below code
        # identically -- except for the random numbers!
        kwargs = self.function_kwargs[self.selected_fingerprinter]
        equivalent_pairs = self.get_equivalent_pairs(population, self.selected_fingerprinter, kwargs)

        if equivalent_pairs and gparameters.mpi.rank == 0:
            if keep_best:
                best = sorted(population, key=lambda individual: individual.fitness)[0].id

            ids = [i.id for i in population]
            equivalent_pairs = [(a.id, b.id) for a, b in equivalent_pairs]
            # disjoint_set_merge will incldue all ids in `ids` as separate entities even if an id is not in any of `equivalent_pairs`
            equivalent_sets = disjoint_set_merge(ids, equivalent_pairs)
            killed = set()
            for equivalent_individuals in equivalent_sets:
                if len(equivalent_individuals) < 2:
                    continue
                if keep_best and best in equivalent_individuals:
                    for x in equivalent_individuals:
                        if x is not best:
                            killed.add(x)
                else:
                    equivalent_individuals = sorted(equivalent_individuals, key=lambda id: population[id].fitness)
                    for x in equivalent_individuals[1:]:
                        killed.add(x)

            check = False
            while len(population) - len(killed) < nkeep:
                rand = random.choice(range(len(killed)))
                rand = MPI.COMM_WORLD.bcast(rand, root=0)
                #killed.pop(rand)
                killed.pop()
                check = True
            if check:  # Make sure each core is killing the same individuals
                all_killed = MPI.COMM_WORLD.allgather(tuple(id for id in killed))
                assert len(set(all_killed)) <= 1

            new_population = [individual for individual in population if individual.id not in killed]
            killed = [population[id] for id in killed]
            population.replace(new_population)
            return killed
        else:
            return []

    @staticmethod
    @parallel
    def get_equivalent_pairs(population, fingerprinter, fingerprinter_kwargs):
        """Returns pairs of individuals that are equivalent.

        Args:
            population (Population): the population
        """
        rank = gparameters.mpi.rank
        ncores = gparameters.mpi.ncores
        pairs_per_core = {r: [] for r in range(ncores)}
        for i, pair in enumerate(combinations(population, 2)):
            pairs_per_core[i % ncores].append(pair)

        equivalent_pairs_by_core = []
        for pair in pairs_per_core[rank]: 
            are_the_same = fingerprinter(*pair, **fingerprinter_kwargs)
            if are_the_same:
                equivalent_pairs_by_core.append(pair)
        if len(equivalent_pairs_by_core) > 0:
            print("Found {} equivalent pairs on rank {}".format(len(equivalent_pairs_by_core), rank))
        count = MPI.COMM_WORLD.allgather(len(equivalent_pairs_by_core))
        all_equivalent_pairs = MPI.COMM_WORLD.allgather(equivalent_pairs_by_core)
        all_equivalent_pairs = [pair for pairs in all_equivalent_pairs for pair in pairs]
        assert sum(count) == len(all_equivalent_pairs)

        return all_equivalent_pairs

    @single_core
    def post_processing(self):
        pass

    @staticmethod
    @functools.wraps(all_close_atom_positions)
    def all_close_atom_positions(individual1, individual2, **kwargs):
        return all_close_atom_positions(individual1, individual2, **kwargs)

    @staticmethod
    @functools.wraps(diversify_module)
    def diversify_module(individual1, individual2, **kwargs):
        return diversify_module(individual1, individual2, **kwargs)

