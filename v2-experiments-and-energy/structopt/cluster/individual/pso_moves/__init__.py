import logging
import functools
from collections import defaultdict

import structopt
from .update_particle import update_particle
from structopt.tools import root, single_core, parallel


class Pso_Moves(object):
    """ """

    @single_core
    def __init__(self, parameters):
        # These variables never change
        self.parameters = parameters
        self.kwargs = self.parameters.update_particles.kwargs


    @single_core
    def move(self, individual, best_swarm, best_particle):
        """ """
        logger = logging.getLogger("default")        
        individual._relaxed = False
        individual._fitted = False
        return self.update_particle(individual, best_swarm, best_particle, self.kwargs.omega, self.kwargs.phi_p, self.kwargs.phi_g)


    @single_core
    def post_processing(self):
        pass


    @staticmethod
    @functools.wraps(update_particle)
    def update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g):
        return update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g)

