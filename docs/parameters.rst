.. _paramters:

Input Parameters
################
The input parameters is a JSON file, which holds a single dictionary. Due to the modular nature of StructOpt, the input file is a dictionary of embedded dictionary, where dictionary key and values often determine the function name and kwargs, respectively. An example Au nanoparticle input file is given below. Each part will be discussed in the following sections.

Example::

  {
    "seed": 0,
    "structure_type": "cluster",
    "generators": {
	"fcc": {"number_of_individuals": 5,
		"kwargs": {"atomlist": [["Au", 55]],
                           "orientation": "100",
			   "cell": [20, 20, 20],
                           "a": 4.08}}
    },
    "fitnesses": {
        "LAMMPS": {"weight": 1.0, 
	           "kwargs": {"use_mpi4py": false,
	                      "MPMD": 0,
	                      "keep_files": true,
	                      "min_style": "cg",
                              "min_modify": "line quadratic",
	                      "minimize": "1e-8 1e-8 5000 10000",
	                      "pair_style": "eam",
	                      "potential_file": "$STRUCTOPT_HOME/potentials/Au_u3.eam",
	                      "thermo_steps": 0,
                              "reference": {"Au": -3.930}}}
    },
    "relaxations": {
        "LAMMPS": {"order": 0,
                   "kwargs": {"use_mpi4py": false,
	                      "MPMD": 0,
	                      "keep_files": true,
	                      "min_style": "cg",
                              "min_modify": "line quadratic",
	                      "minimize": "1e-8 1e-8 5000 10000",
	                      "pair_style": "eam",
	                      "potential_file": "$STRUCTOPT_HOME/potentials/Au_u3.eam",
	                      "thermo_steps": 0}}
    },
    "convergence": {
        "max_generations": 10
    },
    "mutations": {
	"move_surface_atoms": {"probability": 0.0},
        "move_atoms": {"probability": 0.0},
        "move_atoms_group": {"probability": 0.0},
	"rotate_atoms": {"probability": 0.0},
	"rotate_cluster": {"probability": 0.0},
        "rotate_all": {"probability": 0.0},
        "move_surface_defects": {"probability": 1.0}
    },
    "crossovers": {
        "rotate": {"probability": 0.7,
                   "kwargs": {"repair_composition": true}}
    },
    "predators": {
        "fuss": {"probability": 0.9},
        "diversify_module": {"probability": 0.1,
                             "kwargs": {"module": "LAMMPS",
                                        "min_diff": 0.01}}
    },
    "selections": {
        "rank": {"probability": 1.0}
    }
  }

Global Parameters
=================

Global parameters are those that determine how the optimizer should run.

structure_type
++++++++++++++

``structure_type`` (str) is a key parameter for determining operations will be run. Currently, only cluster is supported, but StructOpt is written in a way that *how* the operations are carried out is seperated from the operations themselves. Hence, one can easily write new crossover, mutation, selection, and predator operations that are unique to their structure, test them, and incorporate them inside StructOpt seamlessly.

seed
++++
``seed`` (int): seed for the psuedo-random number generator. Two runs with the exact same input *and* seed should run exactly the same way. Almost all operations use random number generators. See should be an int.

convergence
+++++++++++

``convergence`` (dict): Convergence is a dictionary that determines when to stop the calculation. Currently, the only convergence criteria is ``max_generations``, which is set to an int. For example, the setting below runs the optimizer for 200 generations.

Example::

    "convergence": {
        "max_generations": 200
    }

In the future, more ``convergence`` options will be added.

post_processing
+++++++++++++++

``post_processing`` (dict): Determines what is output as the optimizer is run. Currently, the only option is ``XYZs``, which determines how frequentely the xyz files of each generation should be printed. The rules for this are as follows.

- ``XYZs`` = 0: all generations are kept
- ``XYZs`` > 0: every ``XYZs`` generation is kept
- ``XYZs`` < 0: the last ``XYZs`` generations are kept

The below example is a run where only the last generation is kept. This behavior is by default and encouraged for saving space.
  
Example::

    "post_processing": {
        "XYZs": -1
    }


Generators
==================

Generators are functions for initializing the population. These are pseudo-random generators that depend on the ``seed`` global parameter.

Generators are given as a dictionary entry defined by the ``generators`` key in the input file. The structure of the generators dictionary with *N* desired generators is given below.

Example::

    "generators": {
        generator_1: {"number_of_individuals": n_1,
                      "kwargs": kwargs_1}
        generator_2: {"number_of_individuals": n_2,
                      "kwargs": kwargs_2}
        generator_3: {"number_of_individuals": n_3,
                      "kwargs": kwargs_3}
        ...
        generator_N: {"number_of_individuals": n_N,
                      "kwargs": kwargs_N}
    }

The string for *generator_i*, is the name of the generator one wants to use. The number of individuals that generator should generate is determined by the integer *n_i*. The sum of all *n_i* values determines the total size of the population, which is fixed throughout the run. *kwargs_i* are dictionaries that input the kwargs to the generator function one is using. These will be specific to the function and can be found in their help function, show below.

.. autofunction:: structopt.cluster.individual.generators.ellipsoid

.. autofunction:: structopt.cluster.individual.generators.sphere

.. autofunction:: structopt.cluster.individual.generators.fcc

    
Crossovers
==========

Crossovers are operations for mating two individuals chosen by a selection algorithm. The purpose of the crossover is to intelligently combine (mate) different individuals (parents) in a way to create new individuals (children) that have the features of the current best individuals in the population.

Crossovers are given as a dictionary entry defined by the ``crossovers`` key in the input file. The structure of the crossovers dictionary with *N* desired selections is given below.

Example::

    "crossovers": {
        crossover_1: {"probability": p_1,
                      "kwargs": kwargs_1}
        crossover_2: {"probability": p_2,
                      "kwargs": kwargs_2}
        crossover_3: {"probability": p_3,
                      "kwargs": kwargs_3}
        ...
        crossover_N: {"probability": p_N,
                      "kwargs": kwargs_N}
    }

The string for *crossover_i*,  is the name of the crossover one wants to use. The probability *p_i* is the probability of the crossover occuring if a mate is determined to happen in the population. *p_i* values should sum to 1. *kwargs_i* are dictionaries that input the kwargs to the crossover function one is using. These will be specific to the function and can be found in their help function.

Currently the only crossover in use in the algorithm is the cut-and-splice operator introduced by Deaven and Ho. The description is shown below.

.. autofunction:: structopt.cluster.population.crossovers.rotate

Selections
==========

Selections are operations for choosing which individuals to "mate" in producing new individuals. Individuals are chosen based on their fitness, and different selection functions determine how diverse the children should be.

Selections are given as a dictionary entry defined by the ``selections`` key in the input file. The structure of the selections dictionary with *N* desired selections is given below.

Example::

    "selections": {
        selection_1: {"probability": p_1,
                      "kwargs": kwargs_1}
        selection_2: {"probability": p_2,
                      "kwargs": kwargs_2}
        selection_3: {"probability": p_3,
                      "kwargs": kwargs_3}
        ...
        selection_N: {"probability": p_N,
                      "kwargs": kwargs_N}
    }

The string for *selection_i*,  is the name of the selection one wants to use. The probability *p_i* is the probability of the selection occuring if a mate is determined to happen in the population. *p_i* values should sum to 1. *kwargs_i* are dictionaries that input the kwargs to the selection function one is using. These will be specific to the function and can be found in their help function.

.. autofunction:: structopt.common.population.selections.tournament

.. autofunction:: structopt.common.population.selections.random_selection

.. autofunction:: structopt.common.population.selections.best

.. autofunction:: structopt.common.population.selections.rank

.. autofunction:: structopt.common.population.selections.roulette


Predators
=========

Similar to selections, predators are selection processes that selects individuals based on their fitness. The distinction is that while selections select individuals with positive features to duplicate in children, predators select which individuals to keep in the next generation. Note, this must be done because crossovers and (sometimes) mutations increase the population every generation, and hence each generation requires a predator step.

Predators are given as a dictionary entry defined by the ``predators`` key in the input file. The structure of the predators dictionary with *N* desired predators is given below

Example::

    "predators": {
        predator_1: {"probability": p_1,
                     "kwargs": kwargs_1}
        predator_2: {"probability": p_2,
                     "kwargs": kwargs_2}
        predator_3: {"probability": p_3,
                     "kwargs": kwargs_3}
        ...
        predator_N: {"probability": p_N,
                     "kwargs": kwargs_N}
    }

The string for *predator_i*, is the name of the predator one wants to use. The probability *p_i* is the probability of the predator occuring on every individual in the population. *p_i* values should sum to 1. *kwargs_i* are dictionaries that input the kwargs to the predator function one is using. These will be specific to the function and can be found in their help function.

.. autofunction:: structopt.common.population.predators.best

.. autofunction:: structopt.common.population.predators.roulette

.. autofunction:: structopt.common.population.predators.rank

.. autofunction:: structopt.common.population.predators.fuss

.. autofunction:: structopt.common.population.predators.tournament

.. autofunction:: structopt.common.population.predators.diversify_module
                  
Mutations
=========

Mutations are operations applied to individuals that change its structure and composition. It is a local search operation, though the mutation itself can be written to perform small or larger changes.

Mutations are given as a dictionary entry defined by the ``mutations`` key in the input file. The structure of the mutations dictionary with *N* desired mutations is given below

Example::

    "mutations": {
        "preserve_best": "true" or "false",
        "keep_original": "true" or "false",
        "keep_original_best": "true" or "false,
        mutation_1: {"probability": p_1,
                     "kwargs": kwargs_1}
        mutation_2: {"probability": p_2,
                     "kwargs": kwargs_2}
        mutation_3: {"probability": p_3,
                     "kwargs": kwargs_3}
        ...
        mutation_N: {"probability": p_N,
                     "kwargs": kwargs_N}
    }

The string for *mutation_i*,  is the name of the mutation one wants to use. The probability *p_i* is the probability of the mutation occuring on every individual in the population. *p_i* values should sum to any value between 0 and 1. *kwargs_i* are dictionaries that input the kwargs to the mutation function one is using. These will be specific to the function and can be found in their help function.

In addition to specifying the mutations you want to use, the ``mutations`` dictionary takes three special kwargs: ``preserve_best``, ``keep_original``, and ``keep_original_best``. Setting ``preserve_best`` to ``true``, means the highest fitness individual will **never** be mutated. Setting ``keep_original`` to ``true`` means mutations will be applied to copies of individuals, not the individual itself. This means, the original individual is not changed through a mutation. ``keep_original_best`` applies ``keep_original`` to only the best individual.

The currently implemented mutations are shown below. Note in all functions, the first argument is the atomic structure, which inserted by the optimizer. The user defines all of the other kwargs *after* the first input.

.. autofunction:: structopt.common.individual.mutations.swap_positions

.. autofunction:: structopt.common.individual.mutations.swap_species

.. autofunction:: structopt.common.individual.mutations.rotate_atoms

.. autofunction:: structopt.common.individual.mutations.rotate_all

.. autofunction:: structopt.common.individual.mutations.permutation

.. autofunction:: structopt.common.individual.mutations.rattle

.. autofunction:: structopt.cluster.individual.mutations.move_atoms

.. autofunction:: structopt.cluster.individual.mutations.move_surface_atoms

.. autofunction:: structopt.cluster.individual.mutations.rotate_cluster
                  
.. autofunction:: structopt.cluster.individual.mutations.twist

.. autofunction:: structopt.cluster.individual.mutations.swap_core_shell

.. autofunction:: structopt.cluster.individual.mutations.rich2poor

.. autofunction:: structopt.cluster.individual.mutations.poor2rich

.. autofunction:: structopt.cluster.individual.mutations.move_surface_defects

.. autofunction:: structopt.cluster.individual.mutations.enrich_surface

.. autofunction:: structopt.cluster.individual.mutations.enrich_bulk

.. autofunction:: structopt.cluster.individual.mutations.enrich_surface_defects

.. autofunction:: structopt.cluster.individual.mutations.enrich_surface_facets

Relaxations
===========

Relaxations performs a local relaxation to the atomic structure before evaluating their fitness. This is typically done after crossover and mutation operators are applied.

Relaxations differ than the previous operations in that they require varying amounts of resources. Hence, a subsequent section, Parallelization, will introduce ways to run your job with varying levels of parallel performance.

Relaxations are given as a dictionary entry defined by the ``relaxations`` key in the input file. The structure of these dictionaries is shown below.

Example::

    "relaxations": {
        relaxation_1: {"order": o_1,
                       "kwargs": kwargs_1}
        relaxation_2: {"order": o_2,
                       "kwargs": kwargs_2}
        relaxation_3: {"order": o_3,
                       "kwargs": kwargs_3}
        ...
        relaxation_N: {"order": o_N,
                       "kwargs": kwargs_N}
    }

The string for *relaxation_i*,  is the name of the relaxation one wants to use. The order *o_i* is the order of the relaxation occuring on every individual in the population. *kwargs_i* are dictionaries that input the kwargs to the relaxation function one is using. These will be specific to the function. More details of each relaxation module will be given in the following subsections

LAMMPS
++++++

The LAMMPS relaxation module calls LAMMPS to relax according to some potential. Most of the kwargs can be found from the LAMMPS documentation.

.. autoclass:: structopt.common.individual.relaxations.LAMMPS

The potential files available to use are listed below and are from the default potentials included from LAMMPS. Given a potential, enter in the ``potential_file`` kwarg as "$STRUCTOPT_HOME/potentials/<name>". Note also that different potentials will have different lines of the ``pair_style`` kwarg. If the user would like to use an unavailable potential file, please submit an email to zxu39@wisc.edu, and the potential will be added.

`AlCu.eam.alloy`: Aluminum and copper alloy EAM (Cai and Ye, Phys Rev B, 54, 8398-8410 (1996))

`Au_u3.eam`: Gold EAM (SM Foiles et al, PRB, 33, 7983 (1986))

`ZrCuAl2011.eam.alloy`: Zirconium, copper, and aluminum glass (Howard Sheng at GMU. (hsheng@gmu.edu))

Fitnesses
=========

Fitness evaluates how closely the individual satisfies the minimization criteria. One typical minimization criteria is the stability of a structure, and the formation energy is the fitness. Note, all fitness modules operate so that the *lower* the fitness value the *more* fit it is.

Fitnesses differ than the previous operations in that they require varying amounts of resources. Hence, a subsequent section, Parallelization, will introduce ways to run your job with varying levels of parallel performance.

Fitnesses are given as a dictionary entry defined by ``fitnesses`` key in the input file. The structure of these dictionaries is shown below.

Example::

    "fitnesses": {
        fitness_1: {"weight": w_1,
                     "kwargs": kwargs_1}
        fitness_2: {"weight": w_2,
                     "kwargs": kwargs_2}
        fitness_3: {"weight": w_3,
                     "kwargs": kwargs_3}
        ...
        fitness_N: {"weight": w_N,
                     "kwargs": kwargs_N}
    }


The string for *fitness_i*,  is the name of the fitness one wants to use. The weight *w_i* is the constant to multiply the fitness value returned by the *fitness_i* module. Not that all selections and predators operate on the **total** fitness, which is a sum of each fitness and their weight.  *kwargs_i* are dictionaries that input the kwargs to the fitness function one is using. These will be specific to the function. More details of each fitness module will be given in the following subsections

LAMMPS
++++++

The LAMMPS fitness module calls LAMMPS to calculate the potential energy of the structure. Most of the kwargs can be found from the LAMMPS documentation. In addition, most of the *kwargs* are the same as relaxations, except the fitness module of LAMMPS has a number of normalization options for returning the potential energy. These are described below.

.. autoclass:: structopt.common.individual.fitnesses.LAMMPS

The potential files available to use are listed below and are from the default potentials included from LAMMPS. Given a potential, enter in the ``potential_file`` kwarg as "$STRUCTOPT_HOME/potentials/<name>". Note also that different potentials will have different lines of the ``pair_style`` kwarg. If the user would like to use an unavailable potential file, please submit an email to zxu39@wisc.edu, and the potential will be added.

`AlCu.eam.alloy`: Aluminum and copper alloy EAM (Cai and Ye, Phys Rev B, 54, 8398-8410 (1996))

`Au_u3.eam`: Gold EAM (SM Foiles et al, PRB, 33, 7983 (1986))

`ZrCuAl2011.eam.alloy`: Zirconium, copper, and aluminum glass (Howard Sheng at GMU. (hsheng@gmu.edu))

Parallelization
===============

In addition to the module-specific parameters, each module requires two parallelization entries: ``use_mpi4py`` and ``MPMD_cores_per_structure``. These two entries are mutually exclusive, meaning that only one can be turned on at a time. ``use_mpi4py`` can take two values, ``true`` or ``false`` depending on whether the module should use the `one-structure-per-core <>`_ parallelization.

``MPMD_cores_per_structure`` can be disabled (if ``use_mpi4py`` is ``true``) by setting it to ``0``, but otherwise specifies the number of cores that each process/structure should be allocated within the ``MPI_Comm_spawn_multiple`` command. There are two types of valid values for this parameter: 1) an integer specifying the number of cores per structure, or 2) a string of two integers separated by a dash specifying the minimum and maximum number of cores allowed (e.g. ``"4-16"``). ``MPMD_cores_per_structure`` can also take the value of ``"any"``, and StructOpt will use as many cores as it can to run each individual.

Example::

    "relaxations": {
        "LAMMPS": {"order": 0,
                   "kwargs": {"use_mpi4py": true,
                              "MPMD_cores_per_structure": 0,
                              "keep_files": true,
                              "min_style": "cg\nmin_modify line quadratic",
                              "minimize": "1e-8 1e-8 5000 10000",
                              "pair_style": "eam/alloy",
                              "potential_file": "$STRUCTOPT_HOME/potentials/ZrCuAl2011.eam.alloy",
                              "thermo_steps": 1000}}
    }

