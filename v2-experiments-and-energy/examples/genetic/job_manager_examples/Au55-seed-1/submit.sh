#!/bin/bash

#PBS -N job_manager_examples/Au55-seed-1
#PBS -q morgan2
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -j oe

source activate py35
export PYTHONPATH=$HOME/research/StructOpt_modular/:$PYTHONPATH
export PATH=/share/apps/openmpi-1.10.0_no_ib/bin/:$PATH


cd /home/zxu/research/StructOpt_modular-dev/examples/job_manager_examples/Au55-seed-1

mpirun -n 12 python genetic.py structopt.in.json