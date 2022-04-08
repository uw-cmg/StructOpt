#!/bin/bash
#PBS -N perfect_opt
#PBS -q morgan2
#PBS -l nodes=1:ppn=12,pvmem=1000mb
#PBS -l walltime=96:00:00
##export all environment variables. 
#PBS -V
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
structoptfolder=StructOpt
LD_LIBRARY_PATH=//share/apps/mvapich2/tam_mvapich2-1.9a2/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
/home/zhewen/anaconda2/bin/mpirun python /home/zhewen/$structoptfolder/StructOpt/Optimizer.py structopt_inp.json > out.txt
