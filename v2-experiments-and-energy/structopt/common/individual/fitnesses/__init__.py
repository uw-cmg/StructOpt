from .FEMSIM import FEMSIM
from .LAMMPS import LAMMPS
from .STEM import STEM
from .RDF import RDF
from structopt.tools import root, single_core, parallel


class Fitnesses(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters
        self.modules = []
        self.module_names = []

        for module in self.parameters:
            # Initialize the class that was imported at the top of the file and append it to the modules list
            setattr(self, module, globals()[module](parameters=parameters[module]))
            self.modules.append(getattr(self, module))
            self.module_names.append(module)


    @parallel
    def calculate_fitness(self, individual):
        """Perform the fitness calculations on an individual.

        Args:
            individual (Individual): the individual to evaluate
        """
        fit = 0.0
        # Run each fitness module on the population. Create sorted
        # module list so all cores run modules in the same order
        modules_module_names = [[module, module.__name__.split('.')[-1]] for module in self.modules]
        modules_module_names.sort(key=lambda i: i[1])
        for module, module_name in modules_module_names:
            module_parameters = self.parameters[module_name]
            weight = getattr(module_parameters, 'weight')
            fit += module.calculate_fitness(individual) * weight
        individual._fitness = fit
        self._fitted = True
        return fit


    @single_core
    def post_processing(self):
        pass

