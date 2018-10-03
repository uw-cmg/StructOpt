.. _examples:

Examples
########

The following below are examples of runs you can use to test StructOpt. They exclusively use LAMMPS to relax the structures and calculate its fitness. All of the input files can be found in the StructOpt_modular/examples folder

Running StructOpt
-----------------

StructOpt can be run on a single processor or in parallel. In a single score environment, the command is given below

::
   
   python $STRUCTOPT_HOME/structopt/genetic.py structopt.in.json
   

In a parallel environment with N processors, StructOpt can be run with the following command

::
   
   mpirun -n N python $STRUCTOPT_HOME/structopt/genetic.py structopt.in.json
   

The output will exist in the folder the command was run from

Example 1: cluster/Au55
-----------------------

Example 2: cluster/Au55-parallel
--------------------------------
