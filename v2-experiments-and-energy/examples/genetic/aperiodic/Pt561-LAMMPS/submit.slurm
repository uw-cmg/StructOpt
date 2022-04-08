#!/bin/sh
#SBATCH --partition=stem        # default "univ", if not specified
#SBATCH --time=0-01:30:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1           # require 1 node
#SBATCH --ntasks-per-node=20            # (by default, "ntasks"="cpus")

module load lammps-31Jan14
module load compile/intel
module load mpi/intel/openmpi-1.10.2

python $STRUCTOPT_HOME/structopt/optimizers/genetic.py structopt.in.json 
