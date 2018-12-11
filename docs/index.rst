.. structopt documentation master file, created by
   sphinx-quickstart on Thu May 12 13:35:36 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
   See the following link for documentation syntax:
   https://pythonhosted.org/an_example_pypi_project/sphinx.html

   TODO:  see this link for an example documentation layout:
   http://topik.readthedocs.io/en/latest/api/topik.html


Welcome to StructOpt's documentation!
#####################################

StructOpt is a reverse structure determination toolkit.

What is reverse structure determination?
========================================

Reverse structure determination is a type of structural refinement that iteratively modifies an atomic structure by, for example, moving atoms in the structure with the goal to minimize a function (such as system energy). In atomistic simulations, the positions of the atoms are moved within the model at every step. After the atoms have moved, the structure is evaluated to see how "good" it is. If the structure is "better" than the previous step, the moved atoms are more likely to persist into the next generation. This process is repeated many times until acceptable structure(s) have been generated.

Many different metrics can be used to determine how "good" a structure is, and this is often material-dependent. The average energy per atom is one commonly used metric, and others include fits to experimental data (e.g. S(q) or g(r) data), medium-range order information available via FEM measurements, average coordination number, bond angle constraints, etc.

Different optimization algorithms can be used to minimize the different metrics including Monte Carlo, Genetic Algorithm, and Particle Swarm algorithms.


Overview of StructOpt
=====================

User Documentation
------------------

StructOpt is a structure optimization framework that incorportates multiple forward simulation techniques into its optimization scheme with the goal of identifying stable and realistic atomic structures. It is designed with modularity in mind. Nearly any forward simulation technique that takes an atomic model as input and outputs a fitness value can be integrated into this framework.

This documentation serves as both a user and developer guide for StructOpt. However, parts of this documentation are likely lacking. If you have questions, please post them as an issue on github.

StructOpt serves the purpose of structure refinment for multiple different materials including nanoparticles and metallic glasses and is highly customizable and extendable to new structure types. There are many different types of simulations that can be set up, which requires getting to know the relevent parameters. Examples are included in the github repository and comments via issues on our `github page <https://github.com/uw-cmg/StructOpt>`_ are welcome.

Contents
========

.. toctree::
    :maxdepth: 1

    core_concepts
    setup
    parameters
    outputs
    examples
    parallelism/index
    job_manager
    calculators/index
    why_python
    troubleshooting
    API Reference <structopt/index>

Contributing
============

Bug fixes and error reports are always welcome. We accept PRs and will try to fix issues that have detailed descriptions and are reproducable.

If you have a forward simulation module that you wish to contribute, please make an issue and the correct people will get email notifications so we can respond.


License Agreement
=================

StructOpt is distributed under the `MIT license <https://opensource.org/licenses/MIT>`_, reproduced below:

Copyright (c) 2016 University of Wisconsin-Madison Computational Materials Group

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Index and Search
================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

