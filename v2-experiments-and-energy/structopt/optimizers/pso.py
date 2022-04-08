import os
import sys
import random
import logging

import structopt
import gparameters
from structopt.common.population import Population


class ParticleSwarmOptimization(object):
    """Defines methods to run a particle swarm optimization using the functions in the rest of the library."""

    def __init__(self, population, convergence):
        self.logger = logging.getLogger('default')

        self.population = population
        self.convergence = convergence
        gparameters.generation = 0
        self.converged = False


    def run(self):
        if gparameters.mpi.rank == 0:
            print("Starting main Optimizer loop!")
        while not self.converged:
            self.step()
        if gparameters.mpi.rank == 0:
            print("Finished!")


    def step(self):
        if gparameters.mpi.rank == 0:
            print("Starting generation {}".format(gparameters.generation))
        sys.stdout.flush()
        self._is_best_swarm_updated = False
        if gparameters.generation == 0:
            for id in range(1, len(self.population)):
                self.population[id] = self.population[0].copy()
                self.population[id].id = id
                self.population[id].rattle(stdev=0.5, seed=id)

            self.population.relax()
            fits = self.population.calcualte_fitnesses()
            self.best_swarm = self.population[0].copy()
            self.best_particles = [individual.copy() for individual in self.population]

        self.population.relax()
        fits = self.population.calculate_fitnesses()

        if gparameters.generation > 0:
            for i in range(len(self.population)):
                if fits[i] < self.best_particles[i]._fitness:
                    self.best_particles[i] = self.population[i].copy()
                    if self.best_particles[i]._fitness < self.best_swarm._fitness:
                        self.best_swarm = self.best_particles[i].copy()
                        self._is_best_swarm_updated = True

        updated_population = self.population.run_pso_moves(self.best_swarm, self.best_particles)
        self.population.replace(updated_population)
        self.check_convergence()
        self.post_processing_step()
        gparameters.generation += 1


    def check_convergence(self):
        if gparameters.generation >= self.convergence.max_generations:
            self.converged = True
        else:
            self.converged = False

    def post_processing_step(self):
        if not self._is_best_swarm_updated:
            return
        path = os.path.join(gparameters.logging.path, 'BestXYZs')
        os.makedirs(path, exist_ok=True)
        structopt.io.write_xyz(os.path.join(path, 'generation{}.xyz'.format(gparameters.generation)), self.best_swarm, self.best_swarm._fitness)

    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        pass


if __name__ == "__main__":
    import numpy as np

    parameters = structopt.setup(sys.argv[1])

    random.seed(parameters.seed)
    np.random.seed(parameters.seed)

    population = Population(parameters=parameters)

    with ParticleSwarmOptimization(population=population,
                                   convergence=parameters.convergence
                                   ) as optimizer:
        optimizer.run()

