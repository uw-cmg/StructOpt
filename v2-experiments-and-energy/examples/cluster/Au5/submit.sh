#!/bin/bash

#PBS -N cluster/Au5-with-optimizer
#PBS -q morgan2
#PBS -l nodes=1:ppn=1
#PBS -l walltime=2:00:00
#PBS -j oe

source activate py35
export PYTHONPATH=$HOME/research/StructOpt_modular/:$PYTHONPATH
export PATH=/share/apps/openmpi-1.10.0_no_ib/bin/:$PATH


cd /home/zxu/research/StructOpt_modular/examples/cluster/Au5

mpirun -n 1 python /home/zxu/research/StructOpt_modular/structopt/genetic.py structopt.in.json