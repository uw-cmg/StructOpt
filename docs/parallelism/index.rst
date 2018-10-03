Parallelism
###########

In general, the only parallelized parts of StructOpt are the fitness and relaxation modules that can be plugged in.

StructOpt has two parallelization mechanisms. The first is the simplest case where each structure is assigned to a single core. The core does the significant processing for one structure by processing the module's code. This is optimal when the module does not implement MPI, or the code is relatively fast.

The second parallelization method, called MPMD (via ``MPI_Comm_spawn_multiple``), is a type of advanced dynamic process management but remains relatively easy to use within StructOpt. It allows MPI code to be used within modules and for those modules to be processes on an arbitrary number of cores.

For functions that are only run on the root core (e.g. crossovers and mutations), the `root decorator <https://github.com/uw-cmg/StructOpt_modular/blob/master/structopt/tools/parallel.py>`_ is used on the main ``fitness`` or ``relaxation`` function to broadcast the return value of the function to all cores.



StructOpt acts as a master process ("master program" may be a better word) that runs in Python and uses MPI (via `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_) to communicate between cores. This master process/program makes ``MPI_Comm_spawn_multiple`` calls to C and Fortran programs (which also use MPI). While the C and Fortran processes run, the master python program waits until they are finished. As an example in this section, we will assume StructOpt is using 16 cores to do calculations on 4 structures.

In terms of MPMD parallelization, StructOpt does two primary things:

1. Uses MPI to do preprocessing for the spawning in step (2). ``MPI_Barrier`` is called after this preprocessing to ensure that all ranks have finished their preprocessing before step (2) begins. Note that the preprocessing is distributed across all 16 cores (via the one-core-per structure parallelism using ``mpi4py``), and at the end of the preprocessing the resulting information is passed back to the root rank (e.g. rank == 0).

2. After the preprocessing, the root rank spawns 4 workers, each of which use 4 cores (i.e. all 16 cores are needed to run all 4 processes at the same time). These workers are spawned through either a relaxation or fitness evaluation module, which is done via ``MPI_Comm_spawn_multiple``. These workers can use MPI to communicate within their 4 cores. In the master StructOpt program, only the root rank spawns the C or Fortran subprocesses, and the modules wait until the spawned processes finish before they continue execution.

3. Step (1) and (2) are repeated until the convergence criteria are satisfied.

Cores per Structure Use Cases
-------------------------

* ``ncores == len(population)``: One core per structure

* ``ncores < len(population)``: One core per structure, but all the structure cannot be run at once

* ``ncores > len(population)``: Multiple cores per structures

Unfortunately, it is impossible to predict the number of structures that will be need to be relaxed and fitted after crossovers and mutations have been performed on the population. As a result, all of the above cases are possible (and probable) for any given simulation.


mpi4py: One structure per core
==============================

Main idea:  One structure per core, or multiple structures per core that execute serially in a for-loop. The module must be written in python (or callable from python like LAMMPS through ASE) and implemented directly into StructOpt.

`mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_ allows MPI commands to be run within python. 

Installation
""""""""""""

TODO: Change OpenMPI 1.10.2 to the correct version after the bugfixes have been made. In the meantime, use ``-mca btl tcp,sm,self`` to use TCP rather than infiniband.

Note: mpi4py needs to be installed from source against OpenMPI 1.10.2. Follow the instructions `here <https://media.readthedocs.org/pdf/mpi4py/latest/mpi4py.pdf>`_ under "3.3: Using distutils". In short:

::

    # Setup modules so that `mpi/intel/openmpi` is loaded and `mpiexec` finds that executable
    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-X.Y.tar.gz
    tar -zxf mpi4py-X.Y.tar.gz
    cd mpi4py-X.Y
    python setup.py build
    python setup.py install --user


MPMD: Multiple cores per structure
==================================

Multiple program, multiple data (MPMD) is a form of MPI parallelization where multiple MPI communicators are used synchonously to run multiple MPI processes at the same time. MPMD can be used within ``mpiexec`` by separating each command by colons. Each command is preceded by the ``-n`` option whcih specifies the number of cores to be used for that executable. MPMD can also be used from another MPI master process which calls ``MPI_Comm_spawn_multiple``. This is how StructOpt implements its advanced parallelization techniques to integrate MPI relaxation and fitness programs into its framework. The executable needs to implement MPMD by disconnecting a parent process if it exists (see `here <https://github.com/jjmaldonis/mpi-parallelization/blob/master/spawn_multiple_loop.py>`_ and `here <https://github.com/paul-voyles/femsim-hrmc/blob/master/src/hrmc.f90>`_ for an example parent/child implementation).

