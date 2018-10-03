import logging
import numpy as np

from . import LAMMPS, FEMSIM, STEM
from structopt.tools import root, single_core, parallel
import gparameters


class Fitnesses(object):
    """Holds the parameters for each fitness module and defines a utility function to compute the fitnesses for each fitness module."""

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.modules = [globals()[module] for module in self.parameters]


    @parallel
    def calculate_fitnesses(self, population):
        """Perform the fitness calculations on an entire population.

        Args:
            population (Population): the population to evaluate
        """
        to_fit = [individual for individual in population if not individual._fitted]
        if not to_fit:
            return [individual.fitness for individual in population]

        fitnesses = np.zeros((len(population),), dtype=np.float)
        # Run each fitness module on the population. Create sorted
        # module list so all cores run modules in the same order
        modules_module_names = [[module, module.__name__.split('.')[-1]] for module in self.modules]
        modules_module_names.sort(key=lambda i: i[1])
        for module, module_name in modules_module_names:
            module_parameters = self.parameters[module_name]

            if gparameters.mpi.rank == 0:
                print("Running fitness {} on the entire population".format(module_name))

            fits = module.fitness(population, parameters=module_parameters)

            # Calculate the full objective function with weights
            weight = getattr(module_parameters, 'weight')
            fits = np.multiply(fits, weight)
            fitnesses = np.add(fitnesses, fits)

        # Store the individuals total fitness for each individual and set each individual to
        # unmodified so that the fitnesses won't be recalculated
        for i, individual in enumerate(population):
            individual._fitness = fitnesses[i]
            individual._fitted = True

        self.post_processing(fitnesses)
        return fitnesses

    @single_core
    def post_processing(self, fitnesses):
        logger = logging.getLogger("output")
        logger.info("Total fitnesses for the population: \n{} (rank {})".format(fitnesses, gparameters.mpi.rank))

