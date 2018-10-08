import importlib
import numpy as np
from collections import Counter
from reprlib import recursive_repr as _recursive_repr

import structopt
from ..individual import Individual
from structopt.tools import root, single_core, parallel, allgather
from structopt.tools import SortedDict

POPULATION_MODULES = ['crossovers', 'selections', 'predators', 'fingerprinters', 'fitnesses', 'relaxations', 'mutations', 'pso_moves']

class Population(SortedDict):
    """A list-like class that contains the Individuals and the operations to be run on them."""

    @single_core
    def __init__(self, parameters, individuals=None):
        super().__init__()

        self.parameters = parameters
        self.structure_type = self.parameters.structure_type.lower()
        self.load_modules()

        self._max_individual_id = 0

        if individuals is None:
            # Import the structure type class: e.g from structopt.crystal import Crystal
            # Unfortunately 'from' doesn't seem to work implicitly
            # so a getattr on the module is needed
            module = importlib.import_module('structopt.{}'.format(self.structure_type))
            title = self.structure_type.title()
            if title == "Aperiodic":
                title = "APeriodic"
            Structure = getattr(module, title)

            # Generate/load initial structures
            starting_id = 0
            for generator in sorted(self.parameters.generators.keys()):
                n = self.parameters.generators[generator].number_of_individuals
                for j in range(n):
                    kwargs = self.parameters.generators[generator].kwargs

                    # For the read_xyz, the input is a list of filenames. These need to be
                    # passed as arguments one by one instead of all at once
                    if generator in ['read_xyz', 'read_extxyz']:
                        kwargs = {'filename': kwargs[j]}

                    generator_parameters = {generator: kwargs}

                    structure = Structure(id=starting_id + j,
                                          relaxation_parameters=self.parameters.relaxations,
                                          fitness_parameters=self.parameters.fitnesses,
                                          mutation_parameters=self.parameters.mutations,
                                          pso_moves_parameters=self.parameters.pso_moves,
                                          generator_parameters=generator_parameters)
                    self.add(structure)
                starting_id += n
        else:
            self.update(individuals)

        self.initial_number_of_individuals = len(self)

    def __iter__(self):
        iterable = super().__iter__()
        curr = next(iterable)
        while curr is not StopIteration:
            yield self[curr]
            curr = next(iterable)


    def __reduce__(self):
        'Return state information for pickling'
        inst_dict = vars(self).copy()
        for k in vars(SortedDict()):
            inst_dict.pop(k, None)
        return self.__class__, (), inst_dict or None, None, iter(self.items())


    @_recursive_repr()
    def __repr__(self):
        'od.__repr__() <==> repr(od)'
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))


    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        for name in POPULATION_MODULES:
             if name in state:
                del state[name]
        return state


    def __setstate__(self, other):
        self.__dict__.update(other)
        self.load_modules()


    def load_modules(self):
        importlib.import_module('structopt.{}'.format(self.structure_type))
        for module in self.parameters:
            if module in POPULATION_MODULES and self.parameters[module] is not None:
                Module = importlib.import_module('structopt.{}.population.{}'.format(self.structure_type, module))
                Module = getattr(Module, module.title())(getattr(self.parameters, module))
                setattr(self, module, Module)

    @single_core
    def position(self, individual):
        """Returns the position of the individual in the population."""
        for i, _individual in enumerate(self):
            if _individual is individual:
                return i


    @single_core
    def get_by_position(self, position):
        """Returns the individual at position `position`."""
        for i, individual in enumerate(self):
            if i == position:
                return individual


    @parallel
    def allgather(self, individuals_per_core):
        """Performs an MPI.allgather on self (the population) and updates the
        correct individuals that have been modified based on the inputs from
        individuals_per_core.

        See stuctopt.tools.parallel.allgather for a similar function.
        """
        to_send = [individual for individual in self]

        positions_per_core = {rank: [self.position(individual) for individual in individuals] for rank, individuals in individuals_per_core.items()}

        correct_individuals = allgather(to_send, positions_per_core)
        self.replace(correct_individuals)


    @parallel
    def bcast(self):
        """Performs and MPI.bcast on self."""
        from mpi4py import MPI
        correct_individuals = MPI.COMM_WORLD.bcast([individual for individual in self], root=0)
        self.replace(correct_individuals)


    @single_core
    def add(self, individual):
        """Adds an Individual to the population."""
        assert isinstance(individual, Individual)
        assert individual.id not in [_individual.id for _individual in self]
        self.update([individual])


    @single_core
    def remove(self, individual):
        """Removes `individual` form the population."""
        del self[individual.id]


    @single_core
    def replace(self, a_list):
        """Deletes the current list of individuals and replaces them with the ones in a_list."""
        if self is a_list:
            return None
        self.clear()
        self.update(a_list)


    @single_core
    def update(self, individuals):
        """Overwrites and adds to the population using the `id` attribute of the individuals as a keyword.
        Assigns an id to an individual if it doesn't already have one."""
        for individual in individuals:
            if individual.id is None:
                individual.id = self.get_new_id()
            if individual.id > self._max_individual_id:
                self._max_individual_id = individual.id
        super().update((individual.id, individual) for individual in individuals)


    extend = update


    @single_core
    def get_new_id(self):
        self._max_individual_id += 1
        return self._max_individual_id


    @parallel
    def crossover(self, pairs):
        """Perform crossovers on the population."""
        children = self.crossovers.crossover(pairs)
        return children


    @root
    def mutate(self):
        """Perform mutations on the population."""
        self.mutations.mutate(self)
        return self


    @root
    def run_pso_moves(self, best_swarm, best_particles):
        """Perform PSO moves on the population."""
        self.pso_moves.move(self, best_swarm, best_particles)
        return self


    @parallel
    def calculate_fitnesses(self):
        """Perform the fitness evaluations on the entire population."""

        return self.fitnesses.calculate_fitnesses(self)


    @parallel
    def relax(self):
        """Relax the entire population."""
        self.relaxations.relax(self)


    @parallel
    def apply_fingerprinters(self):
        """Apply fingerprinters on the entire population."""

        self.fingerprinters.select_fingerprinter()
        killed = []
        if self.fingerprinters.selected_fingerprinter is not None:
            killed = self.fingerprinters.remove_duplicates(self, nkeep=self.initial_number_of_individuals, keep_best=self.parameters.fingerprinters.keep_best)
        return killed


    @parallel
    def kill(self):
        """Remove individuals from the population based on a predator scheme."""
        killed = self.__kill()  # Create a new population on the root core
        self.bcast()  # Broadcast the new population
        return killed  # Return the killed individuals on all core

    @root
    def __kill(self):
        """A private method that guarantees the predator is run only on the root but
        allows `self` to be broadcast without returning `self`."""
        self.predators.select_predator()
        killed = self.predators.kill(self, nkeep=self.initial_number_of_individuals)
        return killed

    @root
    def select(self):
        """Select the individuals in the population to perform crossovers on.
        """
        self.selections.select_selection()
        pairs = self.selections.select(self)
        return pairs

