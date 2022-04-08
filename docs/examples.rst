.. _examples:

Examples
########

Examples can be found in the ``examples/`` directory on github.

Running StructOpt
-----------------

StructOpt can be run on a single processor or in parallel using MPI. Depending on the cluster/environment you are using, you may need to load the following modules:

::
   
   module load lammps-31Jan14
   module load compile/intel
   module load mpi/intel/openmpi-1.10.2

StructOpt can be run serially using the following command:

::

   python $STRUCTOPT_HOME/structopt/optimizers/genetic.py structopt.in.json
   

In a parallel environment with *N* processors, StructOpt can be run with the following command:

::
   
   mpirun -n N python $STRUCTOPT_HOME/structopt/optimizers/genetic.py structopt.in.json
   

The output will exist in the folder the command was run from.
