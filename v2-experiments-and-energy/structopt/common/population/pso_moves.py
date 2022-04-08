from structopt.tools import root, single_core, parallel


class Pso_Moves(object):
    """ """

    @single_core
    def __init__(self, parameters):
        self.parameters = parameters


    @single_core
    def move(self, population, best_swarm, best_particles):
        """Update the entire population.
        Args:
            population (Population): the population to update
        """
        for i, individual in enumerate(population):
            individual.pso_moves.move(individual, best_swarm, best_particles[i])
        return population

    @single_core
    def post_processing(self):
        pass

