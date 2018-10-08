import os
import sys
import logging
import time

import structopt
import gparameters
import structopt.utilities
import structopt.postprocessing
from structopt.common.population import Population
from structopt.tools.convert_time import convert_time


class GeneticAlgorithm(object):
    """Defines methods to run a genetic algorithm optimization using the functions in the rest of the library."""

    def __init__(self, population, convergence, post_processing):
        self.logger = logging.getLogger('default')

        self.population = population
        self.convergence = convergence
        self.post_processing = post_processing

        gparameters.generation = 0
        self.converged = False

        self.timing = {'step': [],
                       'fitness': [],
                       'relax': [],
                       'selection': [],
                       'crossover': [],
                       'mutation': [],
                       'predator': []}

    def run(self):
        if gparameters.mpi.rank == 0:
            print("Starting main Opimizer loop!")
        while not self.converged:
            self.step()
        if gparameters.mpi.rank == 0:
            print("Finished running GA!")


    def step(self):
        t_step_0 = time.time()

        if gparameters.mpi.rank == 0:
            print('')
            print("Starting generation {}".format(gparameters.generation))
        sys.stdout.flush()
        if gparameters.generation > 0:
            fits = [individual._fitness for individual in self.population]

            t_selection_0 = time.time()
            parents = self.population.select(fits)
            self.timing['selection'].append(time.time() - t_selection_0)

            t_crossover_0 = time.time()
            children = self.population.crossover(parents)
            self.population.extend(children)
            self.timing['crossover'].append(time.time() - t_crossover_0)

            t_mutation_0 = time.time()
            mutated_population = self.population.mutate()
            self.population.replace(mutated_population)
            self.timing['mutation'].append(time.time() - t_mutation_0)
        else:
            self.timing['selection'].append(0)
            self.timing['crossover'].append(0)
            self.timing['mutation'].append(0)

        t_relax_0 = time.time()
        self.population.relax()
        self.timing['relax'].append(time.time() - t_relax_0)

        t_fitness_0 = time.time()
        fits = self.population.fitness()
        self.timing['fitness'].append(time.time() - t_fitness_0)

        t_predator_0 = time.time()
        self.population.kill(fits)
        self.timing['predator'].append(time.time() - t_predator_0)

        self.timing['step'].append(time.time() - t_step_0)

        self.check_convergence()
        if gparameters.mpi.rank == 0:
            self.post_processing_step()
        gparameters.generation += 1

    def check_convergence(self):
        if gparameters.generation >= self.convergence.max_generations:
            self.converged = True
        else:
            self.converged = False

    def post_processing_step(self):
        # Save the fitnesses for each individual
        fitness_logger = logging.getLogger('fitness')
        for individual in self.population:
            line = 'Generation {}, Individual {}:'.format(gparameters.generation, individual.id)
            for module in individual.fits:
                line += ' {}: {}'.format(module, individual.fits[module])
            fitness_logger.info(line)

        # Save the XYZ file for each individual.
        for individual in self.population:
            path = os.path.join(gparameters.logging.path, 'XYZs')
            os.makedirs(path, exist_ok=True)
            path = os.path.join(path, 'generation{}'.format(gparameters.generation))
            os.makedirs(path, exist_ok=True)
            individual.write(os.path.join(path, 'individual{}.xyz'.format(individual.id)))
        structopt.postprocessing.clear_XYZs(self.post_processing.XYZs, gparameters.generation, gparameters.logging.path)

        # Save the genealogy
        tags = ['' for _ in self.population]
        for i, individual in enumerate(self.population):
            tags[i] = '{ctag}{id}{mtag}'.format(ctag=individual.crossover_tag or '', id=individual.id, mtag=individual.mutation_tag or '')
            individual.crossover_tag = None
            individual.mutation_tag = None
        genealogy_logger = logging.getLogger('genealogy')
        genealogy_logger.info(' '.join(tags))

        # Save the times
        timing_logger = logging.getLogger('timing')
        timing_logger.info('')
        timing_logger.info('Generation {} (cumulative) timing information'.format(gparameters.generation))
        for operation in ['selection', 'crossover', 'mutation',
                          'relax', 'fitness', 'predator', 'step']:
            t, t_unit = convert_time(self.timing[operation][-1])
            t_cum, t_cum_unit = convert_time(sum(self.timing[operation]))
            timing_logger.info('{:10s}: {:4.2f} {} ({:4.2f} {})'.format(operation, t, t_unit, t_cum, t_cum_unit))

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass


if __name__ == "__main__":
    import structopt
    import random
    import numpy as np

    parameters = structopt.setup(sys.argv[1])
    random.seed(parameters.seed)
    np.random.seed(parameters.seed)

    population = Population(parameters=parameters)

    with GeneticAlgorithm(population=population,
                          convergence=parameters.convergence,
                          post_processing=parameters.post_processing) as optimizer:
        optimizer.run()

