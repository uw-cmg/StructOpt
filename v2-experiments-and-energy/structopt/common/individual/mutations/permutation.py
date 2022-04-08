import random

def permutation(individual):
    """Swaps the chemical symbol between two elements

    Parameters
    ----------
    individual : Individual
        An individual or atoms object.
    """

    a1 = individual[random.randint(0, individual.get_number_of_atoms() - 1)]
    opts = [i for i in individual if i.symbol != a1.symbol]
    a2 = opts[random.randint(0,len(opts) - 1)]
    a1.symbol, a2.symbol = a2.symbol, a1.symbol

    return
