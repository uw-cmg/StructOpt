from .LAMMPS import LAMMPS
from .STEM import STEM
from .hard_sphere_cutoff import hard_sphere_cutoff
from structopt.tools import root, single_core, parallel


class Relaxations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.modules = []

        for module in self.parameters:
            # Initialize the class that was imported at the top of the file and append it to the modules list
            parameters = self.parameters[module]
            setattr(self, module, globals()[module](parameters=parameters))
            self.modules.append(getattr(self, module))


    @parallel
    def relax(self, individual, generation=None):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """
        for module in self.modules:
            module.relax(individual, generation)
        individual._relaxed = True
        individual._fitted = False
        return None


    @single_core
    def post_processing(self):
        pass

