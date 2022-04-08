import sys

from core import DataExplorer


dir = sys.argv[1]

de = DataExplorer(dir)
last_generation = len(de) - 1

#print(de.parameters)

ID = 81
GENERATION = 5

# Access via generation-population-id
indiv = de[GENERATION][ID]
print(id(indiv))
print(indiv)
print(indiv.LAMMPS)
print(indiv.STEM)
print(indiv.fitness)

# Access history via id
ih = de.get_historical_individual(ID)
print(ih.created_on)
print(ih.killed_on)
for id_ in ih:
    ind = ih[id_]
    print(id_, ind, ind.crossover_tag, ind.mutation_tag, ind.get_parents(de))

histories = set()
for pop in de:
    for id_, ind in pop.items():
        histories.add(ind.history)

ages = set()
for ih in histories:
    age = ih.killed_on - ih.created_on
    ages.add((age, ih))
ages = sorted(ages, key=lambda t: t[0], reverse=True)
print(ages[0])

#print("History:")
#indiv.print_history(de.generations)

#crossovers, mutations = de.generations[-1][id_].count_modifiers(de.generations)
#print("Crossovers:", crossovers)
#print("Mutations:", mutations)

indiv.load_structure()
print(indiv.get_positions()[:5])
print(indiv.relaxations)
