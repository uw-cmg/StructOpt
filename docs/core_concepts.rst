.. _core_concetps:

Core Concepts
#############

Detailed information about StructOpt can be found in our paper: [in the submission process]

Overview and General Workflow
-----------------------------

StructOpt uses a Genetic Algorithm to optimize a set of atomic structures according to a customizable objective function (aka cost function).

Genetic Algorithm
=================

A genetic algorithm utilizes a population of structures rather than a single individual. A genetic algorithm, or evolutionary algorithm, is conceptually similar to genetic Darwinism where animals are replaced by "individuals". In StructOpt an "individual" is an atomic model. A population of atomic models is first generated. Given this population, pairs of individuals are mated (aka crossed over) by selecting different volumes of different models and joining them. Crossovers always produce two children, one for each section of the models combined together. The offspring are added to the population. After the mating scheme has finished, single individuals can "mutate" (i.e. moving atoms in a unique way) to add new "genes" to the population's gene pool. After the atoms have been moved via crossovers and mutations, the structures are relaxed. Finally, each structure is run though a series of "fitness" evaluations to determine how "fit to survive" it is, and the population is then reduced to its original size based on a number of optional selection criteria. This process is repeated many times.

In summary:

1. Generate initial structures
2. Locally relax structures
3. Calculate fitness values (e.g. energies) of each structure
4. Remove some individuals from the population based on their fitness value
5. Perform crossovers and selected individuals to generate offspring for the next generation
6. Perform mutations on the selected individuals in the current population and offspring for the next generation
7. Repeat steps 2-6 until the convergence criteria are met

The relaxation and fitness evaluations are only performed on individuals that have been modified via crossovers and mutations. This avoids recomputing these expensive calculations for individuals that were unchanged during the generation's crossover/mutation scheme.

During crossovers, the offspring are collected into a list. After all crossovers have been completed, these offspring are added to the entire population. Each individual in the entire population then has a chance to be mutated. There will therefore be a variable number of modified individuals that will need to be relaxed and fit during each generation. The number of modified individuals can only be predicted by using the probability of mutation and crossover.
