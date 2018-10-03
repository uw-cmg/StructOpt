import logging

from . import LAMMPS
from . import STEM
from . import hard_sphere_cutoff
from structopt.tools import root, single_core, parallel
import gparameters


class Relaxations(object):
    """Holds the parameters for each relaxation module and defines a utility function to run the relaxations for each relaxation module."""

    @single_core
    def __init__(self, parameters):

        # Store the relaxation modules in a list based on their specified order
        self.parameters = parameters
        modules = [globals()[module] for module in self.parameters]
        orders = [self.parameters[module]['order'] for module in self.parameters]
        modules_orders = list(zip(modules, orders))
        modules_orders = sorted(modules_orders, key=lambda modules_orders: modules_orders[1])
        try:
            self.modules = list(zip(*modules_orders))[0]
        except IndexError:
            self.modules = []

    @parallel
    def relax(self, population):
        """Relax the entire population using all the input relaxation methods.

        Args:
            population (Population): the population to relax
        """
        logger = logging.getLogger("default")
        to_relax = [individual for individual in population if not individual._relaxed]
        logger.info("Found {} individuals to relax on core {}: {}".format(len(to_relax), gparameters.mpi.rank, to_relax))
        if not to_relax:
            return

        for i, module in enumerate(self.modules):
            if gparameters.mpi.rank == 0:
                print("Running relaxation {} on the entire population".format(module.__name__.split('.')[-1]))
            parameters = self.parameters[module.__name__.split('.')[-1]]
            module.relax(population, parameters=parameters)

        for individual in population:
            individual._relaxed = True

        return


    @single_core
    def post_processing(self):
        pass

