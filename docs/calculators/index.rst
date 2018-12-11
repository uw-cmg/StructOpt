Relaxation and Fitness Modules
##############################

LAMMPS
------

Installation
============

Follow the `standard installation instructions <http://lammps.sandia.gov/doc/Section_start.html>`_.

Create an environment variable called ``LAMMPS_COMMAND`` that points to the serial LAMMPS executable after installation.

`Package Documentation <http://lammps.sandia.gov/>`_


VASP
----

VASP currently cannot be run within StructOpt. We have done quite a bit of work to get it running, but it isn't working yet. If you'd like to work on this, please post an issue on github.


FEMSIM
------

Installation
============

Fork and clone the `repository <https://github.com/paul-voyles/femsim-hrmc>`_ from github.

Using OpenMPI 1.10.2 compilers, follow the instructions to compile femsim.

Create an environment variable called ``FEMSIM_COMMAND`` pointing to the newly created ``femsim`` executable.


`Package Documentation <https://github.com/paul-voyles/femsim-hrmc>`_


STEM
----

The STEM software is included in the StructOpt repository.

References:  http://pubs.acs.org/doi/abs/10.1021/acsnano.5b05722


Creating Your Own Module
------------------------

Any forward simulation that takes an atomic model as input and outputs a "fitness" value that can be interpreted as a measure of "goodness" of a structure can be integrated into StructOpt. Contact the developers by making an issue on github to get in touch with us.

