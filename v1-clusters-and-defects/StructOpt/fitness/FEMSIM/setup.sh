#!/bin/bash

git clone https://github.com/paul-voyles/femsim-hrmc.git
cd femsim-hrmc
# If this following command fails make sure the mpif90 being used (on bardeen) is:
# FC = /share/apps/mvapich2/1.2.1-p1/gcc_ifort/bin/mpif90
make femsim -C src
