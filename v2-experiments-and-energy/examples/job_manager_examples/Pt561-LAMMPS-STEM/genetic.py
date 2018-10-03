import sys, os
import logging
import shutil

import structopt
import structopt.utilities
from structopt.common.population import Population

class GeneticAlgorithm(object):
    """Defines methods to run a genetic algorithm optimization using the functions in the rest of the library."""

    def __init__(self, population, convergence, post_processing, adaptation):
        self.logger = logging.getLogger('default')

        self.population = population
        self.convergence = convergence
        self.post_processing = post_processing
        self.adaptation = adaptation

        self.generation = 0
        self.converged = False


    def run(self):
        if logging.parameters.rank == 0:
            print("Starting main Opimizer loop!")
        while not self.converged:
            self.step()
        if logging.parameters.rank == 0:
            print("Finished running GA!")


    def step(self):
        if logging.parameters.rank == 0:
            print('')
            print("Starting generation {}".format(self.generation))
        sys.stdout.flush()
        if self.generation > 0:
            fits = [individual._fitness for individual in self.population]
            parents = self.population.select(fits)
            children = self.population.crossover(parents)
            self.population.extend(children)
            mutated_population = self.population.mutate()
            self.population.replace(mutated_population)

        self.population.relax()
        fits = self.population.fitness()
        self.population.kill(fits)
        self.check_convergence()
        if logging.parameters.rank == 0:
            self.post_processing_step()
        structopt.utilities.adapt(self.adaptation, self.population, self.generation)
        logging.parameters.generation += 1
        self.generation += 1


    def check_convergence(self):
        if self.generation >= self.convergence.max_generations:
            self.converged = True
        else:
            self.converged = False

    def post_processing_step(self):
        # Save the fitnesses for each individual
        fitness_logger = logging.getLogger('fitness')
        for individual in self.population:
            line = 'Generation {}, Individual {}:'.format(self.generation, individual.id)
            for module in individual.fits:
                line += ' {}: {}'.format(module, individual.fits[module])
            fitness_logger.info(line)

        # Save the XYZ file for each individual.
        for individual in self.population:
            path = os.path.join(logging.parameters.path, 'XYZs')
            os.makedirs(path, exist_ok=True)
            path = os.path.join(path, 'generation{}'.format(self.generation))
            os.makedirs(path, exist_ok=True)
            individual.write(os.path.join(path, 'individual{}.xyz'.format(individual.id)))

        structopt.utilities.clear_XYZs(self.post_processing.XYZs, self.generation, logging.parameters.path)

        # Save the genealogy
        tags = ['' for _ in self.population]
        for i, individual in enumerate(self.population):
            tags[i] = '{ctag}{id}{mtag}'.format(ctag=individual.crossover_tag or '', id=individual.id, mtag=individual.mutation_tag or '')
            individual.crossover_tag = None
            individual.mutation_tag = None
        genealogy_logger = logging.getLogger('genealogy')
        genealogy_logger.info(' '.join(tags))

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


if __name__ == "__main__":
    import sys
    import structopt
    import random
    import numpy as np
    from genetic import GeneticAlgorithm

    parameters = structopt.setup(sys.argv[1])
    random.seed(parameters.seed)
    np.random.seed(parameters.seed)
    
    population = Population(parameters=parameters)

    with GeneticAlgorithm(population=population,
                          convergence=parameters.convergence,
                          post_processing=parameters.post_processing,
                          adaptation=parameters.adaptation) as optimizer:    
        optimizer.run()

