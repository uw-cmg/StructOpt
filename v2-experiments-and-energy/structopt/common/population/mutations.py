import numpy as np

from structopt.tools import root, single_core, parallel


class Mutations(object):
    """ """

    @single_core
    def __init__(self, parameters):
        """
        Keywords in parameters:
        keep_original : Add mutated individuals to the population without removing the individual they were mutated from.
        keep_original_best : If the individual with the lowest fitness is mutated, do not remove its un-mutated parent from the population.
        """
        self.parameters = parameters
        self.keep_original = parameters.get('keep_original', False)
        self.keep_original_best = parameters.get('keep_original_best', False)

    @single_core
    def mutate(self, population):
        if self.keep_original or self.keep_original_best:
            fits = {individual.id: (individual.fitness if individual._fitted else np.inf) for individual in population}
            min_fit_id = min(fits, key=fits.get)

        # Make sure not to edit the population in the for loop!
        to_remove = []
        to_add = []
        for individual in population:
            individual.mutations.select_mutation()

            if individual.mutations.selected_mutation is not None:
                # Duplicate the individual and reset some values
                mutated = individual.copy()
                mutated.mutated_from = individual.id
                mutated.mutations.selected_mutation = individual.mutations.selected_mutation
                individual.mutations.selected_mutation = None

                # Perform the mutation
                mutated.mutate(select_new=False)

                # Replace the individual with the mutated one
                if not self.keep_original and not (self.keep_original_best and individual.id == min_fit_id):
                    to_remove.append(individual)
                to_add.append(mutated)

        for individual in to_remove:
            population.remove(individual)
        for mutated in to_add:
            population.add(mutated)

        return population


    @single_core
    def post_processing(self):
        pass

