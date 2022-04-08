#PBS -S /bin/bash
#PBS -m be
#PBS -q cpu16mem64
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=72:00:00
#PBS -N JOB_NAME

# available queues: cpu16mem64, cpu20mem128
# replace the 1 in select= with the number of nodes (but use only 1 - no working IB)
# replace the 16 in ncpus and mpiprocs with 20 depending on the queue
# specify walltime in the HH:MM:SS format
# replace JOB_NAME with your desired jobe name

cd $PBS_O_WORKDIR 

#module load vasp
module load intel
module load impi
module load python3

# use prun with no argumets, no matter what MPI environment is used
python $STRUCTOPT_HOME/structopt/optimizers/genetic.py Al_497.in.json > out 
