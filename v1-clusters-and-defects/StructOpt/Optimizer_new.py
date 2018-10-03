import importlib
import numpy as np

class EmptyStructure(object):
    pass

class Optimizer(object):
    def __init__(self, modules, weights, relaxation_method):
        self.step_number = 0

        self.args = self.read_inputs()
        for k, v in self.args.items():
            setattr(self, k, v)

        self.modules = []
        self.weights = []

        for i, module in enumerate(modules):
            imported = importlib.import_module("StructOpt.fitness.{module}.{module}_eval".format(module=module))
            cls = getattr(imported, '{cls_name}_eval'.format(cls_name=module))  # Get class from module
            self.add_module(cls(), weights[i])

        self.relaxation_method = relaxation_method

        # Define GA parameters
        self.individuals = self.initialize_structures()
        # self... = ...


    def read_inputs(self):
        return {}


    def initialize_structures(self, restart=False):
        # Create structures (or reload if restart==True)
        return [EmptyStructure(), EmptyStructure(), EmptyStructure()]


    def objective_function(self):
        fitnesses = np.array([0.0 for _ in self.individuals])
        for module, weight in zip(self.modules, self.weights):
            fitness = np.array(module.evaluate_fitness(self.individuals))
            fitnesses += weight*fitness
        return fitnesses


    def move_atoms(self):
        """Modifies all the individuals in self.individuals"""
        return self.individuals


    def relax_structures(self):
        """Relaxes all the individuals in self.individuals"""
        return self.individuals


    def add_module(self, module, weight):
        self.modules.append(module)
        self.weights.append(weight)


    def step(self):
        self.strucutures = self.move_atoms()  # Move atoms
        self.strucutures = self.relax_structures()
        self.step_number += 1

        for module in self.modules:
            module.update_parameters()

        fitnesses = self.objective_function()

        # Optional
        for i, individual in enumerate(self.individuals):
            individual.fitness = fitnesses[i]

        return fitnesses


    def run(self):
        self.initialize_structures()
        while not self.converged:
            if self.step_number >= 100:
                break
            self.step()


    @property
    def converged(self):
        """Returns whether or not the simulation has converged"""
        return False
