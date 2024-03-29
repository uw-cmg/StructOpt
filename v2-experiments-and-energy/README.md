[![Documentation Status](https://readthedocs.org/projects/structopt/badge/?version=latest)](http://structopt.readthedocs.org/en/latest/)

StructOpt is a Structure Optimization suite designed for materials with complicated structure. It can incorportate the results from any number of experimental data, provided that experimental data can be simulated. StructOpt's goal is to identify atomic structures that are energetically stable and agree with experimental data. It is designed with modularity in mind, and encourages simplicity in both its codebase and usage without sacrificing powerful functionality. Nearly any forward simulation technique that takes an atomic model as input and outputs a fitness value can be integrated into its optimizer.

StructOpt's most developed optimizer is the genetic algorithm, but particle swarm and monte carlo functionality are available. The user can also develop their own optimizers that simply use the material-modification API to modify and transform structructres. In this way, a new optimizer can take advantage of the complex transformations that are already in use by the genetic algorithm, for example.

StructOpt serves the purpose of structure refinment for multiple different materials including nanoparticles, defects, and metallic glasses. As such, it is highly customizable and extendable, and there are many different types of simulations that can be set up.

### Citing StructOpt

The manuscript is in the submission process.

A snapshot of the code at the time of publication can be found here: https://github.com/uw-cmg/StructOpt/releases/tag/0.1
