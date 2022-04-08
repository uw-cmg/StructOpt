import logging
import random
import ase
from importlib import import_module
import numpy as np

import structopt
from structopt.tools import root, single_core, parallel
from .generate_velocities.random_velocities import random_velocities

class Individual(ase.Atoms):
    """An abstract base class for a structure."""

    @single_core
    def __init__(self, id=None,
                 load_modules=True,
                 relaxation_parameters=None, fitness_parameters=None,
                 mutation_parameters=None,
                 pso_moves_parameters=None,
                 generator_parameters=None, **kwargs):

        self.id = id
        self.mutation_tag = None
        self.crossover_tag = None
        self.relaxation_parameters = relaxation_parameters
        self.fitness_parameters = fitness_parameters
        self.mutation_parameters = mutation_parameters
        self.generator_parameters = generator_parameters
        self.pso_moves_parameters = pso_moves_parameters
        self._fitted = False
        self._relaxed = False
        self._fitness = None
        self._Q_l = np.array([])

        cls_name = self.__class__.__name__.lower()
        if load_modules:
            self.load_modules()
            self.has_modules = True
        else:
            self.has_modules = False

        # Initialize the ase.Atoms structure
        super().__init__()

        # Due to the complexity to most of the generator functions, most of the
        # generators will be classes, not functions.  Making them classes
        # will also allow for class inheritance to eliminate potential duplicate
        # code when using similar functions for generating clusters, bulk, surfaces, etc.
        # However, read_xyz is a function and must be called uniquely.

        if self.generator_parameters is not None:
            self.generate()

        else:
            return None
        random_velocities(self)


    def __eq__(self, other):
        try:
            return self._fitness == other._fitness
        except AttributeError:
            if self._fitness is None and other._fitness is None:
                return True
            else:
                return False


    def __ne__(self, other):
        try:
            return self._fitness != other._fitness
        except AttributeError:
            if self._fitness is None and other._fitness is None:
                return False
            else:
                return True


    def __lt__(self, other):
        try:
            return self._fitness < other._fitness
        except AttributeError:
            if self._fitness is None:
                return False
            elif other._fitness is None:
                return True


    def __le__(self, other):
        try:
            return self._fitness <= other._fitness
        except AttributeError:
            if self._fitness is None and other._fitness is None:
                return True
            elif self._fitness is None:
                return False
            elif other._fitness is None:
                return True


    def __gt__(self, other):
        try:
            return self._fitness > other._fitness
        except AttributeError:
            if self._fitness is None:
                return True
            elif other._fitness is None:
                return False


    def __ge__(self, other):
        try:
            return self._fitness >= other._fitness
        except AttributeError:
            if self._fitness is None and other._fitness is None:
                return True
            elif self._fitness is None:
                return True
            elif other._fitness is None:
                return False


    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries. The unpickled object WILL NOT have these attributes at all!
        for name in ['fitnesses', 'relaxations', 'mutations', 'pso_moves']:
            if name in state:
                del state[name]
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.__dict__.update(state)
        if self.has_modules:
            self.load_modules()


    def __str__(self):
        return '<Individual {}>'.format(self.id)
    __repr__ = __str__


    @parallel
    def load_modules(self):
        """Loads the relevant modules."""
        cls_name = self.__class__.__name__.lower()

        # Load in the appropriate functionality
        if self.fitness_parameters is not None:
            fitnesses = import_module('structopt.{}.individual.fitnesses'.format(cls_name))
            self.fitnesses = fitnesses.Fitnesses(parameters=self.fitness_parameters)
            for name in self.fitnesses.module_names:
                if not hasattr(self, name):
                    setattr(self, name, None)
        else:
            self.fitnesses = None

        if self.mutation_parameters is not None:
            mutations = import_module('structopt.{}.individual.mutations'.format(cls_name))
            self.mutations = mutations.Mutations(parameters=self.mutation_parameters)
        else:
            self.mutations = None

        if self.relaxation_parameters is not None:
            relaxations = import_module('structopt.{}.individual.relaxations'.format(cls_name))
            self.relaxations = relaxations.Relaxations(parameters=self.relaxation_parameters)
        else:
            self.relaxations = None

        if self.pso_moves_parameters is not None:
            pso_moves = import_module('structopt.{}.individual.pso_moves'.format(cls_name))
            self.pso_moves = pso_moves.Pso_Moves(parameters=self.pso_moves_parameters)
        else:
            self.pos_moves = None


    @property
    @single_core
    def velocities(self):
        return self.get_velocities()


    @property
    @single_core
    def positions(self):
        return self.get_positions()

    def wrap(self, center=(0.5, 0.5, 0.5), pbc=None, eps=1e-7):
        if pbc is None:
            pbc = self.pbc
        self.set_positions(ase.geometry.wrap_positions(
            self.positions, self.cell, pbc, center, eps))

    @property
    @single_core
    def fits(self):
        return {module: getattr(self, module) for module in self.fitnesses.module_names}


    @parallel
    def mutate(self, select_new=True):
        """Mutate an individual.

        Args:
            individual (Individual): the individual to mutate
        """
        if select_new:
            self.mutations.select_mutation()
        self.mutations.mutate(self)


    @parallel
    def relax(self):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        self.relaxations.relax(self)


    @parallel
    def calculate_fitness(self):
        """Perform the fitness calculations on an individual.

        Args:
            individual (Individual): the individual to evaluate
        """
        return self.fitnesses.calculate_fitness(self)


    @property
    def fitness(self):
        """The total fitness of the individual."""
        return self._fitness


    @single_core
    def generate(self):
        """Generate an individual using generator_kwargs parameter. By defualt
        it extends the current atoms object"""

        cls_name = self.__class__.__name__.lower()

        assert len(self.generator_parameters) == 1
        generator_name = list(self.generator_parameters.keys())[0]
        generator_module = import_module('structopt.{}.individual.generators.{}'.format(cls_name, generator_name))
        if not hasattr(generator_module, generator_name):
            generator_module = import_module('structopt.common.individual.generators.{}'.format(generator_name))
        generator = getattr(generator_module, generator_name)
        kwargs = self.generator_parameters[generator_name]
        atoms = generator(**kwargs)
        self.extend(atoms)
        self.set_cell(atoms.get_cell())
        self.set_pbc(atoms.get_pbc())

        return


    @single_core
    def get_atom_indices_within_distance_of_atom(self, atom_index, distance):
        dists = self.get_distances(atom_index, slice(None, None, None))
        return np.where(dists < distance)


    @single_core
    def get_nearest_atom_indices(self, atom_index, count):
        dists = self.get_distances(atom_index, slice(None, None, None))[0]
        return np.argsort(dists)[:count]


    @single_core
    def copy(self, include_atoms=True):
        """Return a copy."""
        new = self.__class__(id=None,
                             load_modules=True,
                             relaxation_parameters=self.relaxation_parameters,
                             fitness_parameters=self.fitness_parameters,
                             mutation_parameters=self.mutation_parameters,
                             pso_moves_parameters=self.pso_moves_parameters,
                             generator_parameters=None)
        if include_atoms:
            new.arrays = {}
            for name, a in self.arrays.items():
                new.arrays[name] = a.copy()
        else:
            new.clear()
        new.set_cell(self.get_cell())
        new.set_pbc(self.get_pbc())
        new.mutation_tag = self.mutation_tag
        new.crossover_tag = self.crossover_tag
        new._fitted = self._fitted
        new._relaxed = self._relaxed
        new._fitness = self._fitness
        new._Q_l = self._Q_l
        if self.fitnesses is not None:
            for module_name in self.fitnesses.module_names:
                setattr(new, module_name, getattr(self, module_name, None))
        return new

    @single_core
    def clear(self):
        del self[:]
