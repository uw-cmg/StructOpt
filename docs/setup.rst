Installation and Setup
######################

StructOpt is written in Python 3 and as such requires a working Python 3 installation. We recommend setting up an Anaconda virtual environment exclusively for StructOpt.


Python Libraries
----------------

::

    conda install numpy
    conda install scipy
    pip install ase
    pip install natsorted
    # Install mpi4py from source (below)

mpi4py
======

On Madison's ACI cluster:

::

    module load compile/intel
    module load mpi/intel/openmpi-1.10.2

Follow `these <http://mpi4py.readthedocs.io/en/stable/install.html#using-distutils>`_ instructions:  

::

    wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-X.Y.tar.gz
    tar -zxf mpi4py-X.Y.tar.gz
    cd mpi4py-X.Y
    python setup.py build
    python setup.py install --user

You can test your installation by following `these <http://mpi4py.readthedocs.io/en/stable/install.html#testing>`_ instructions.


Installing StructOpt
--------------------

To get the code, fork and clone the StructOpt repository or download the zip `here <https://github.com/uw-cmg/StructOpt_modular>`_. Add the location of the StructOpt folder (e.g. ``$HOME/repos/StructOpt``) to your ``PATH`` environment variable.

Create an environment variable called ``STRUCTOPT_HOME`` with the same folder location as you added to your path.


Additional Modules
------------------

Depending on the type of calculation you wish to run, you will need to install specific relaxation and fitness modules. VASP and LAMMPS are two examples that function as both a relaxation and fitness. They can both move atoms and relax the structure, and they output an energetic term that quantifies how "good" or fit the input structure is. These modules should be installed using their standard installation procedures unless specified otherwise in the `Modules <>`_ section.
