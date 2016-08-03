from structopt.common.individual import Individual


class Cluster(Individual):
    """ A stucture containing a non-periodic cluster of atoms. Relaxation algorithms such as LAMMPS and VASP require
        a periodic structure, so when run though these types of algorithms, the cluster is embedded in a larger box.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.set_pbc([False, False, False])

