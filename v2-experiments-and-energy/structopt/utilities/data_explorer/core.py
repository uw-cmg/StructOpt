import os
import sys
import json
import re
import importlib
import weakref
from collections import defaultdict, Counter
import numpy as np
import warnings

from .common import lazy, lazyproperty

from structopt.io import read_xyz
from structopt.tools.dictionaryobject import DictionaryObject


def parse_gene(gene):
    gene_tester = re.compile('(?P<id>\d+)?(?P<cross>c\D+\(\d+\+\d+\))?(?P<mut>m\D+\(\d+\))?')
    matches = gene_tester.search(gene).groups()
    id, crossover, mutation = matches
    return int(id), crossover or None, mutation or None


class DataExplorer(object):
    def __init__(self, dir):
        self.genealogy_file = os.path.join(dir, 'genealogy.log')
        self.fitnesses_file = os.path.join(dir, 'fitnesses.log')
        self.output_file = os.path.join(dir, 'output.log')
        self._icache = {}
        self._created_killed = None

    @lazyproperty
    def generations(self):
        return Generations(self.genealogy_file, self)

    @lazyproperty
    def parameters(self):
        parameters = []
        start_flag = False
        for line in open(self.output_file):
            if 'Current parameters:' in line:
                start_flag = True
                continue
            if start_flag:
                if 'INFO : {' in line:
                    line = line.split('INFO : ')[1]
                parameters.append(line)
            if line == '}\n':
                start_flag = False
                break
        parameters = ''.join(parameters)
        parameters = DictionaryObject(json.loads(parameters))
        sys.modules['gparameters'] = parameters
        return parameters

    @lazyproperty
    def _fitnesses(self):
        fitnesses = {}
        for line in open(self.fitnesses_file):
            firsthalf = re.compile(".+ : INFO : Generation (\d+), Individual (\d+):(.*)")
            secondhalf = re.compile("([\w]+): ([-\d]+.\d+|inf)")
            generation, id, fitness_str = firsthalf.findall(line)[0]
            fitnesses[int(id)] = {fit: float(value) for fit, value in secondhalf.findall(fitness_str)}
        return fitnesses

    def _get_fitness(self, id):
        fits = self._fitnesses[id]
        weights = {name: module.weight for name, module in self.parameters.fitnesses.items()}
        fitness = sum(fits[module]*weights[module] for module in fits)
        return fitness

    def _get_module_fitness(self, id, module):
        return self._fitnesses[id][module]

    def __getitem__(self, index):
        return self.generations[index]

    def __len__(self):
        return len(self.generations)

    def get_individual(self, id):
        """Returns the Individual object for the given id or the id if not found."""
        created, killed = self._get_created_killed(id)
        if created == np.inf or created == -np.inf:
            return id
        return self.generations[created][id]

    def _get_created_killed(self, id):
        """
        Returns
        -------
            tuple<int, int> : (created_on, killed_on)
        """
        if self._created_killed is None:
            created_killed = defaultdict(lambda: [np.inf, -np.inf])
            for generation, population in enumerate(self):
                for id_ in population.keys():
                    if generation < created_killed[id_][0]:
                        created_killed[id_][0] = generation
                    if generation > created_killed[id_][1]:
                        created_killed[id_][1] = generation + 1
            self._created_killed = created_killed
        return self._created_killed[id]


class Generations(object):
    def __init__(self, genealogy_filename, dataexplorer):
        self._dataexplorer = weakref.ref(dataexplorer)
        generations = open(genealogy_filename).readlines()
        self._data = [line.strip().split(' : INFO : ')[-1] for line in generations]

    def __getitem__(self, index):
        if isinstance(self._data[index], Population):
            return self._data[index]
        else:
            generation, data = self._data[index].split(": ")  # Split at : in "Generation n: ..."
            generation = int(generation.split(" ")[-1])  # Get n from "Generation n"
            assert generation == index
            self._data[index] = Population(generation, data, self._dataexplorer())
            return self._data[index]

    def __len__(self):
        return len(self._data)


class Population(dict):
    def __init__(self, generation, individuals_from_genealogy_file, dataexplorer):
        self._generation = generation
        self._dataexplorer = weakref.ref(dataexplorer)
        self._data = individuals_from_genealogy_file
        self._loaded = False

    @property
    def generation(self):
        return self._generation

    def __hash__(self):
        return self.generation

    def _load(self):
        self._data = self._data.split()
        self._data = [parse_gene(ind) for ind in self._data]
        self._data = {id: (c, m) for id, c, m in self._data}
        super().__init__({id: None for id in self._data if id not in self})
        self._loaded = True

    def __getitem__(self, key):
        if key not in self or super().__getitem__(key) is None:
            if not self._loaded:
                self._load()
                if key not in self._data:
                    raise KeyError(key)
            ct, mt = self._data[key]
            self[key] = Individual(key, self._dataexplorer(), ct, mt)
        return super().__getitem__(key)

    def __len__(self):
        return len(self._data)

    def keys(self):
        if not self._loaded:
            self._load()
        return super().keys()

    def items(self):
        if not self._loaded:
            self._load()
        for id in self:
            self[id]  # Access the id to call __getitem__ and make sure the individual is loaded
        return super().items()

    def values(self):
        if not self._loaded:
            self._load()
        for key in self:
            self[key]
        return super().values()

    def __iter__(self):
        if not self._loaded:
            self._load()
        return super().__iter__()


class Individual(object):
    def __init__(self, id, dataexplorer, crossover_tag=None, mutation_tag=None):
        self._id = id
        self._dataexplorer = weakref.ref(dataexplorer)
        self._created_on = None
        self._killed_on = None
        self._loaded = False
        self._structure = None

        self.mutation_tag = mutation_tag
        self.crossover_tag = crossover_tag
        if self.mutation_tag is not None:
            parent = re.compile('m(\w+)\((\d+)\)')
            tag, parent = parent.search(self.mutation_tag).groups()
            self.mutation_tag = tag
            self._mutated_from = int(parent)
        if self.crossover_tag is not None:
            parents = re.compile('c(\w+)\((\d+)\+(\d+)\)')
            tag, parent1, parent2 = parents.search(self.crossover_tag).groups()
            self._parent1 = int(parent1)
            self._parent2 = int(parent2)
            self.crossover_tag = tag

    def __hash__(self):
        return id(self)

    @property
    def id(self):
        return self._id

    @property
    def generations(self):
        return list(range(self.created_on, self.killed_on))

    @property
    def parents(self):
        ps = {}
        if self.mutation_tag is not None:
            ps["mutation"] = self._dataexplorer().get_individual(self._mutated_from)
        if self.crossover_tag is not None:
            ps["crossover"] = [self._dataexplorer().get_individual(self._parent1), self._dataexplorer().get_individual(self._parent2)]
        return ps

    @property
    def created_on(self):
        if self._created_on is not None:
            return self._created_on
        c, k = self._dataexplorer()._get_created_killed(self.id)
        self._created_on = c
        self._killed_on = k
        return self._created_on

    @property
    def killed_on(self):
        if self._killed_on is not None:
            return self._killed_on
        c, k = self._dataexplorer()._get_created_killed(self.id)
        self._created_on = c
        self._killed_on = k
        return self._killed_on

    @lazyproperty
    def fitness(self):
        return self._dataexplorer()._get_fitness(self.id)

    @lazyproperty
    def LAMMPS(self):
        return self._dataexplorer()._get_module_fitness(self.id, "LAMMPS")

    @lazyproperty
    def STEM(self):
        return self._dataexplorer()._get_module_fitness(self.id, "STEM")

    @lazyproperty
    def FEMSIM(self):
        return self._dataexplorer()._get_module_fitness(self.id, "FEMSIM")

    def __repr__(self):
        tag = '{id}{ctag}{mtag}'.format(
                                         id=self.id,
                                         ctag="c{t}({p1}+{p2})".format(t=self.crossover_tag, p1=self._parent1, p2=self._parent2) if self.crossover_tag else '',
                                         mtag="m{t}({p})".format(t=self.mutation_tag, p=self._mutated_from) if self.mutation_tag else '',
                                        )
        return "<Individual {}>".format(tag)
        #return "<Individual {} @{}>".format(tag, self.generations)

    __str__ = __repr__

    def load_structure(self, filename=None):
        parameters = self._dataexplorer().parameters
        if "generators" in parameters:
            parameters.pop("generators")
        self.structure_type = parameters.structure_type.lower()
        module = importlib.import_module('structopt.{}'.format(self.structure_type))
        Structure = getattr(module, self.structure_type.title())
        if filename is None:
            filename = os.path.join(self._dataexplorer().parameters.logging.path, 'modelfiles', 'individual{}.xyz'.format(self.id))
        generator_parameters = {"read_xyz": {"filename": filename}}
        self._structure = Structure(id=self.id,
                              relaxation_parameters=parameters.relaxations,
                              fitness_parameters=parameters.fitnesses,
                              mutation_parameters=parameters.mutations,
                              pso_moves_parameters=parameters.pso_moves,
                              generator_parameters=generator_parameters)
        self._loaded = True
        return self._structure

    def __getattr__(self, key):
        if not self._loaded:
            self.load_structure()
            if hasattr(self._structure, key):
                return getattr(self._structure, key)
        return super().__getattribute__(key)

    def __len__(self):
        return len(self._structure)

    @lazyproperty
    def ancestry(self):
        """
        Returns
        -------
            list<Individual> : All ancestors of this individual.
        """
        ancestry = []
        parents = []
        if "mutation" in self.parents:
            parents.append(self.parents["mutation"])
        if "crossover" in self.parents:
            for p in self.parents["crossover"]:
                parents.append(p)
        ancestry.extend(parents)
        for parent in parents:
            if parent is not None and not isinstance(parent, int):
                ancestry.extend(parent.ancestry)

        if any(isinstance(individual, int) for individual in ancestry):
            warnings.warn("WARNING: At least one parent individual was killed in the same generation it was created, and we therefore do not have data for it. The ancestry is there truncated.")

        return ancestry

    def get_modifier_counts(self):
        ancestry = [individual for individual in self.ancestry if not isinstance(individual, int)]
        counter = Counter()
        for individual in ancestry:
            counter[individual.mutation_tag] += 1
            counter[individual.crossover_tag] += 1
        if None in counter:
            del counter[None]
        return counter

