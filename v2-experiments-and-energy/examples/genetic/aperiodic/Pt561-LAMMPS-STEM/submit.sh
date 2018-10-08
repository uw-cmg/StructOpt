#!/bin/bash

#PBS -N structopt/LAMMPS+STEM/Pt561-cuboctahedron/run23/s0
#PBS -q morgan2
#PBS -l nodes=1:ppn=12
#PBS -l walltime=18:00:00
#PBS -j oe

source activate py35
export PYTHONPATH=$HOME/research/StructOpt_modular/:$PYTHONPATH
export PATH=/share/apps/openmpi-1.10.0_no_ib/bin/:$PATH


cd /home/zxu/research/PtMo-nanoparticles/structopt/LAMMPS+STEM/Pt561-cuboctahedron/run23/s0

mpirun -n 12 python genetic.py structopt.in.json