#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH -t 1-00:00:00
#SBATCH --partition=pre

export SLURM_NNODES; export SLURM_NTASKS=`echo "$SLURM_NNODES*20" | bc`

module load compile/intel
module load mpi/intel/openmpi-1.10.2


date

#/usr/mpi/intel/openmpi-1.10.2/bin/mpirun -np $SLURM_NTASKS python $STRUCTOPT_HOME/structopt/optimizers/genetic.py  structopt.in.json
#mpirun.openmpi -np $SLURM_NTASKS python $STRUCTOPT_HOME/structopt/optimizers/genetic.py  structopt.in.json
python $STRUCTOPT_HOME/structopt/optimizers/genetic.py TiO2_600.in.json
#mpirun -n 2 $STRUCTOPT_HOME/structopt/optimizers/genetic.py TiO2_600.in.json
  
date
