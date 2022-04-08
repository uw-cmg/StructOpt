from ase.neighborlist import NeighborList

from structopt.tools import root, single_core, parallel
import gparameters


class hard_sphere_cutoff(object):
    """A relaxation module to ensure atoms in an individual are not too close together.
    This is often a preliminary relaxation before LAMMPS for VASP to ensure the models do not explode.
    """

    @single_core
    def __init__(self, parameters, cutoff=1.0):
        # These variables never change
        self.cutoff = cutoff
        print(cutoff)
        self.parameters = parameters


    @single_core
    def relax(self, individual):
        """Relaxes the individual using a hard-sphere cutoff method.
        Args:
            individual (Individual):  the individual to relax
        """
        rank = gparameters.mpi.rank
        print("Relaxing individual {} on rank {} with hard-sphere cutoff method".format(individual.id, rank))
        elem = individual.get_chemical_symbols()
        radii = []
        for e in elem:
            if e == 'O':
                radii.append(1.4)  # ionic radii of O2-
            elif e == 'Ti':
                radii.append(0.605) # ionic radii of Ti4+

        nl = NeighborList(radii, bothways=True, self_interaction=False)
        nl.update(individual)

        ntries = 0
        modified = True
        while modified and ntries < 100:
            modified = False
            for atom in individual:
                indices, offsets = nl.get_neighbors(atom.index)
                if 'Ti' in atom.symbol and 'Ti' in individual[indices].symbols:
                    cutoff = 2.0
                elif 'O' in atom.symbol and 'O' in individual[indices].symbols:
                    cutoff = 2.0
                else:
                    cutoff = 2.0
                for neigh in indices:
                    if individual.get_distance(atom.index, neigh) < cutoff:
                        individual.set_distance(atom.index, neigh, cutoff, fix=0.1)
                        modified = True
            nl.update(individual)
            individual.wrap()
            ntries += 1
        if ntries == 100:
            print("WARNING! Iterated through the hard-sphere cutoff relaxation 100 times and it still did not converge!")

