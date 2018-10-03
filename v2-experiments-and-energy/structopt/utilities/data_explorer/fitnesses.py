

def load_fitnesses(filename):
    generations = open(filename).readlines()
    generations = [line.strip().split(' : INFO : ')[-1] for line in generations]
    generations = [line.split() for line in generations]
    return fitnesses

