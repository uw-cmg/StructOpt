[![Documentation Status](https://readthedocs.org/projects/structopt/badge/?version=latest)](http://structopt.readthedocs.org/en/latest/)

StructOpt
=========

StructOpt is a Structure Optimization suite designed to identify atomic structures that are energetically stable and consistent with experimental data. Manuscripts published before 2018 referencing StructOpt refer to `v1-clusters-and-defects`. Manuscripts published in or after 2018 refer to `v2-experiments-and-energy`. See the subdirectories for more information.

ReadTheDocs currently only supports documentation in the head of the repository, not in subdirectories, so the `docs` directory contains the documentation for `v2`.

StructOpt is written in Python 3 and as such requires a working Python 3 installation. We recommend setting up an Anaconda virtual environment exclusively for StructOpt.

## TOC
- [Install](#Install)
- [LAMMPS](#LAMMPS)
- [Run](#Run-code)

## Install
### 0.  Create Anaconda virtual environment with python 3.6 
Install Python Libraries 
```bash
conda install numpy
conda install scipy
conda install diffpy.srreal (https://github.com/diffpy/diffpy.srreal)
pip install pymatgen
pip install ase
pip install natsorted
```
### 1. Check openmpi is properly installed by running  ```mpirun --version```   
```bash
sudo apt install libopenmpi-dev
sudo apt-get install -y openmpi-bin
```
### 2. install mpi4py from source against the openmpi environment
Install mpi4py by running command ```pip install --no-cache-dir mpi4py```.    
If run into error when using mpi4py, check if the mpi4py is installed against the mpi version correctly. See the post https://stackoverflow.com/questions/55129738/centos-7-undefined-symbol-ompi-mpi-logical8 .  
Or install mip4py from source: the mpi4py packages are avaliable at: https://bitbucket.org/mpi4py/mpi4py/downloads/

```bash
download
    curl -0 https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
    or:
    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
install
    tar -zxf mpi4py-3.0.0.tar.gz
    cd mpi4py-3.0.0
    python setup.py build
    python setup.py install --user
    https://mpi4py.readthedocs.io/en/stable/install.html#using-distutils
```
### 3. Install lammps and compile aganist the openmpi environment
```bash
Download from https://github.com/lammps/lammps/tree/stable_3Mar2020  

install  

unzip lammps-stable_3Mar2020.zip
cd lammps-stable_3Mar2020/src
make ps       #check the current installed packages
make yes-user-meamc yes-manybody yes-user-reaxc
make mpi
make serial
```

### 4. Install intel parallel studio  (omit it if you got openmpi on the cluster)
Run following script to install  

```bash
tar xvzf l_compxe_2015.3.187.tgz
cd l_compxe_2015.3.187
./install.sh
```

By default, the installation path is $HOME/intel. During the install step, you need to enter the path of the license file.  
After installation, source the environmental variables. Insert the following lines to your ~/.bashrc.  

```bash
# intel parallel studio
source $HOME/intel/bin/compilervars.sh intel64

# link mpifort to use ifort instead of gfortran
# https://stackoverflow.com/questions/32335519/gfortran-is-called-instead-of-mpif90
export I_MPI_F90=ifort
```
###  5. Install FEMSIM_HRMC
Download from https://github.com/paul-voyles/femsim-hrmc  
You can install femsim following instruction from their github repo by following script  

```bash
make hrmc –C src
make femsim –C src
```

I have some trouble in following their instruction because my mpif90 do not point to a valid gfortran binary. So I have rewrite the hrmc.c and makefile to build serial version. Use the [modified codes](https://github.com/mengjun930/FEMSIM-HRMC) and run  

```bash
make hrmc serial -C src
make femsim serial -C src
```

### 6. Set environmental variables for structopt  

Set the these path according to your install configuration and insert them to ~/.bashrc  

```bash
export LAMMPS_COMMAND=/PATH/TO/lmp
export FEMSIM_COMMAND=/PATH/TO/femsim
export PYTHONPATH=/PATH/TO/StructOpt/v2-experiments-and-energy
export STRUCTOPT_HOME=/PATH/TO/StructOpt/v2-experiments-and-energy
```

## Run code

### Set number of cores for lammps & femsim job  
Several jobs (lammps or femsim) are run in parallel using python's process pool (multiprocessing.Pool). The maximum number of jobs is controlled by relaxations.LAMMPS.njobs and fitnesses.FEMSIM.njobs in [input json file](input/anatase221_0408/TiO2_48.in.json), respectively.  
Each femsim job is run by a single core.  
Each lammps could run on multiple cores using mpirun. The number of cores can be set by ```ncores``` variable in v2-experiments-and-energy/structopt/common/crossmodule/lammps.py.  

### Then run  

```bash
python v2-experiments-and-energy/structopt/optimizers/genetic.py input/alpha0_input/TiO2_600.in.json
```

## Notes
### Checklist for lammps relaxation 
- atom type should be in the same order between structure file & & mass & pair_coeff.  

- Difference between inputs for COMB and MEAM potentail  
    - atom_style charge (for comb) vs atomic (for meam)  
    - while MEAM potential do not need to specify the atomic mass because it can read them from its potential file. The COMB potential needs users to specify the atomic mass explicitly.  
    - comb needs to perform qeq by adding ```fix fix_qeq all qeq/comb 1 0.01```
- Minimization with COMB potential
    - The convergence and outcoming structure of minimization using comb potential depends on the initial charge assigned to each atom. By performing some trial calculations, I found it's better to assign initial charge of -0.5 and 1.0 to O and Ti, respectively. Assigning charge (0, 0) leads to failed minization on some trail calculations. While assigning (-1.5, 3.0), which is closer to final charge distribution (-1.6, 3.2), would give slightly higher energy comparing with (-0.5, 1.0).  
    - relax minimization criterion ```minimize 0.0 1e-4 10000 10000```  
    - relax qeq step ```fix fix_qeq all qeq/comb 10 0.01``` 


