import numpy as np
from scipy.interpolate import interp1d
import os

from ase.calculators.lammpsrun import LAMMPS as lammps

from structopt.tools import root, single_core, parallel
from structopt.tools.dictionaryobject import DictionaryObject
import gparameters


class LAMMPS(object):
    """LAMMPS class for running LAMMPS on a single individual. Takes
    a dictionary, where the key: value are the parameters for running LAMMPs.

    Parameters
    ----------
    min_style : str
        The minimization scheme for running LAMMPS. See LAMMPS doc.
    min_modify : str
        Parameters for min_style energy minimization algorithm.
        See LAMMPS doc.
    minimize : str
        Convergence criteria for minimization algorithm. Note for 
        fitness values, the last two values are set to 0, so no
        relaxation is done. See LAMMPS doc.
    pair_style : str
        Type of potential used. See LAMMPS doc.
    potential_file : str
        The path to the potential_file. Should be absolute.
    thermo_steps : int
        How much output to print of thermodynamic information.
        If set to 0, only the last step is printed.See LAMMPS doc.
    keep_file : bool
        Will keep all of the LAMMPS input and output files for each
        individual. Use with caution.
    reference : dict
        Reference energies of the particle. These are values to subtract
        from the values returned by LAMMPS. Given as a dictionary of
        {sym : E} pairs, where sym is a str denoating the
        the element, while E is the value to be subtracted per sym. This is
        typically the pure component formation energy calculated with LAMMPS.
        Note since this is merely a fixed subtraction, should not change the
        performance in constant composition runs.
    """


    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        self.output_dir = gparameters.logging.path


    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @single_core
    def calculate_fitness(self, individual):
        # Don't rerun lammps if:
        # 1) the individual is unmodified
        # 2) the energy has already been calculated via the relaxation
        if individual._relaxed and 'LAMMPS' in individual.relaxations.parameters and individual.LAMMPS is not None:
            E = individual.LAMMPS
        else:
            print("Individual {} did not have a value for .LAMMPS or it was modified".format(individual.id))
            calcdir = os.path.join(self.output_dir, 'fitness/LAMMPS/generation{}/individual{}'.format(gparameters.generation, individual.id))
            rank = gparameters.mpi.rank

            calc = lammps(self.parameters.kwargs, calcdir=calcdir)
            individual.set_calculator(calc)
            try:
                # We will manually run the lammps calculator's calculate.
                #  Normally calc.calculate would get run with default arguments via:
                #  ase.get_potential_energy -> lammps.get_potential_energy -> lammps.update -> lammps.calculate
                #  but we want to run it with a custom trajectory file output location, so we manually call calculate.
                #  Then, when ase calls calculate, it won't run because it's already been finished.
                trj_file = os.path.join(gparameters.logging.path, "modelfiles", "individual{}.trj".format(individual.id))
                calc.calculate(individual, trj_file=trj_file)
                E = individual.get_potential_energy()
                print("Finished calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, rank))
            except RuntimeError:
                E = np.inf
                print("Error calculating fitness of individual {} on rank {} with LAMMPS".format(individual.id, rank))

        E = self.reference(E, individual)
        E = self.normalize(E, individual)
        individual.LAMMPS = E
        return E


    @single_core
    def reference(self, E, individual):
        """References the energy of the cluster to a reference energy"""
        if 'reference' not in self.parameters:
            return E

        E_refs = self.parameters.reference

        if type(self.parameters.reference) is float:
            return E - self.parameters.reference

        elif isinstance(self.parameters.reference, dict):
            symbols = individual.get_chemical_symbols()
            E_ref = [E_refs[symbol] for symbol in symbols]
            return E - sum(E_ref)

        else:
            raise TypeError('LAMMPS parameter "reference" must be either float or dict')

    @single_core
    def normalize(self, E, individual):
        if 'normalize' not in self.parameters:
            return E

        norms = self.parameters.normalize

        if 'natoms' in norms and norms['natoms'] is True:
            E /= len(individual)

        if 'particle' in norms:
            n = len(individual)
            f = np.poly1d(norms['particle'])
            E -= f(np.log(n))

        if 'alloy' in norms:
            element = norms['alloy']['element']
            xs = norms['alloy']['x']
            Es = norms['alloy']['E']

            linear = np.poly1d(np.polyfit([xs[0], xs[-1]], [Es[0], Es[-1]], 1))
            Es = np.array(Es) - linear(xs)

            f = interp1d(xs, Es)

            x = individual.get_chemical_symbols().count(element) / len(individual)
            E -= f(x)

        return E

