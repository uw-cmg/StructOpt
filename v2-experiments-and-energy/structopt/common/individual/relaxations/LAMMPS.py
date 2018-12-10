import os
import numpy as np

from structopt.common.crossmodule.lammps import LAMMPS as lammps
from structopt.tools import root, single_core, parallel
from structopt.aperiodic.individual.mutations.move_surface_atoms import move_surface_atoms
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
        Convergence criteria for minimization algorithm. See LAMMPS doc.
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
    repair : bool
        Determines whether to run an algorithm to make sure no atoms
        are in "space". Atoms can be in space due to a mutation or
        crossover that results in a large force that shoots the atom
        outside of the particle.
    """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters.kwargs
        if hasattr(gparameters, 'logging'):
            self.output_dir = gparameters.logging.path
        else:
            self.output_dir = '.'

    @single_core
    def get_command(self, individual):
        raise NotImplementedError


    @parallel
    def relax(self, individual):
        """Relax an individual.

        Args:
            individual (Individual): the individual to relax
        """

        calcdir = os.path.join(self.output_dir, 'relaxation/LAMMPS/generation{}/individual{}'.format(gparameters.generation, individual.id))
        rank = gparameters.mpi.rank
        print("Relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        calc = lammps(self.parameters, calcdir=calcdir)
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
            print("Finished relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))
        except RuntimeError:
            E = np.inf
            print("Error relaxing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        individual.LAMMPS = E

        if 'repair' in self.parameters and self.parameters['repair']:
            E = self.repair(individual, gparameters.generation)
            if E is not None:
                individual.LAMMPS = E

        return

    @parallel
    def repair(self, individual, generation):
        """Repairs an individual. Currently takes isolated atoms moves them next to
        a non-isolated atom"""

        rank = gparameters.mpi.rank

        n_moves = move_surface_atoms(individual, max_natoms=1.0, move_CN=3)
        if n_moves == 0:
            return

        if generation is not None:
            calcdir = os.path.join(os.getcwd(), '{}/relaxation/LAMMPS/generation{}/individual{}-repair')
            calcdir = calcdir.format(self.output_dir, generation, individual.id)
        else:
            calcdir = None

        calc = lammps(self.parameters, calcdir=calcdir)
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
            print("Finished repairing individual {} on rank {} with LAMMPS".format(individual.id, rank))
        except RuntimeError:
            E = np.inf
            print("Error repairing individual {} on rank {} with LAMMPS".format(individual.id, rank))

        return E

