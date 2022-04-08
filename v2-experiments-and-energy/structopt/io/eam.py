# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2014) James Kermode, King's College London
#                  Lars Pastewka, Karlsruhe Institute of Technology
#                  Adrien Gola, Karlsruhe Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# ======================================================================


from collections import namedtuple

import numpy as np
try:
    from scipy import interpolate
except:
    print('Warning: No scipy')
    interpolate = False

import os

###

EAMParameters = namedtuple('EAMParameters', 'symbols atomic_numbers '
                           'atomic_masses lattice_constants crystal_structures '
                           'number_of_density_grid_points '
                           'number_of_distance_grid_points '
                           'density_grid_spacing distance_grid_spacing '
                           'cutoff')

###

def read_eam(eam_file, kind="eam/alloy"):
    """
    Read an eam alloy lammps format file and return the tabulated data and parameters
    http://lammps.sandia.gov/doc/pair_eam.html
    
    Parameters
    ----------
      eam_file : string
                      eam alloy file name 
      kind : string
             kind of EAM file to read (supported eam,eam/alloy,eam/fs)
    Returns
    -------
      source : string
          Source informations or comment line for the file header
      parameters : list of tuples
                [0] - array of str - atoms (ONLY FOR eam/alloy and eam/fs, EMPTY for eam)
                [1] - array of int - atomic numbers
                [2] - array of float -atomic masses
                [3] - array of float - equilibrium lattice parameter
                [4] - array of str - crystal structure
                [5] - int - number of data point for embedded function
                [6] - int - number of data point for density and pair functions
                [7] - float - step size for the embedded function
                [8] - float - step size for the density and pair functions
                [9] - float - cutoff of the potentials
      F : array_like
          contain the tabulated values of the embedded functions
          shape = (nb atoms, nb of data points)
      f : array_like
          contain the tabulated values of the density functions
          shape = (nb atoms, nb of data points)
      rep : array_like
          contain the tabulated values of pair potential
          shape = (nb atoms,nb atoms, nb of data points)
    """

    eam = open(eam_file, 'r').readlines()
    if kind == "eam":
        # reading first comment line as source for eam potential data
        source = eam[0].strip()
        # -- Parameters -- #
        atnumber = int(eam[1].split()[0])
        atmass = float(eam[1].split()[1])
        crystallatt = float(eam[1].split()[2])
        crystal = eam[1].split()[3]
        Nrho = int(eam[2].split()[0])       # Nrho (number of values for the embedding function F(rho))
        Nr = int(eam[2].split()[2])         # Nr (number of values for the effective charge function Z(r) and density function rho(r))
        drho = float(eam[2].split()[1])     # spacing in density space
        dr = float(eam[2].split()[3]) # spacing in distance space
        cutoff = float(eam[2].split()[4])
        parameters = EAMParameters(np.zeros(1),atnumber, atmass,crystallatt,crystal, Nrho,Nr, drho, dr, cutoff)
        # -- Tabulated data -- #
        data = np.loadtxt(eam_file, dtype="float", skiprows = 3).flatten()
        F = data[0:Nrho]
        f = data[Nrho:Nrho+Nr]
        rep = data[Nrho+Nr:2*Nr+Nrho]
        return source,parameters, F,f,rep
    elif kind == "eam/alloy":
        # reading 3 first comment lines as source for eam potential data
        source = eam[0].strip()+eam[1].strip()+eam[2].strip()
        # -- Parameters -- #
        atoms = eam[3].strip().split()[1:]
        nb_atoms = len(atoms)
        Nrho = int(eam[4].split()[0])       # Nrho (number of values for the embedding function F(rho))
        Nr = int(eam[4].split()[2])         # Nr (number of values for the effective charge function Z(r) and density function rho(r))
        drho = float(eam[4].split()[1])     # spacing in density space
        dr = float(eam[4].split()[3]) # spacing in distance space
        cutoff = float(eam[4].split()[4])
        atnumber,atmass,crystallatt,crystal = np.empty(nb_atoms,dtype=int),np.empty(nb_atoms),np.empty(nb_atoms),np.empty(nb_atoms).astype(np.str)
        for i in range(nb_atoms):
            # Fixme: The following lines assume that data occurs in blocks of
            # homogeneous width. This can break.
            l = len(eam[6].strip().split())
            row = int(5+i*((Nr+Nrho)/l+1))
            atnumber[i] = int(eam[row].split()[0])
            atmass[i] = float(eam[row].split()[1])
            crystallatt[i] = float(eam[row].split()[2])
            crystal[i] = str(eam[row].split()[3])
        parameters = EAMParameters(atoms,atnumber,atmass,crystallatt,crystal,Nrho,Nr,drho,dr,cutoff)
        # -- Tabulated data -- #
        F,f,rep,data = np.empty((nb_atoms,Nrho)),np.empty((nb_atoms,Nr)),np.empty((nb_atoms,nb_atoms,Nr)),np.empty(())
        eam = open(eam_file,'r')
        [eam.readline() for i in range(5)]
        for i in range(nb_atoms):
            eam.readline()
            data = np.append(data,np.fromfile(eam,count=Nrho+Nr, sep=' '))
        data = np.append(data,np.fromfile(eam,count=-1, sep=' '))
        data = data[1:]
        for i in range(nb_atoms):
            F[i,:] = data[i*(Nrho+Nr):Nrho+i*(Nrho+Nr)]
            f[i,:] = data[Nrho+i*(Nrho+Nr):Nrho+Nr+i*(Nrho+Nr)]
        interaction = 0
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                if j < i :
                    rep[i,j,:] = data[nb_atoms*(Nrho+Nr)+interaction*Nr:nb_atoms*(Nrho+Nr)+interaction*Nr+Nr]
                    rep[j,i,:] = data[nb_atoms*(Nrho+Nr)+interaction*Nr:nb_atoms*(Nrho+Nr)+interaction*Nr+Nr]
                    interaction+=1
            rep[i,i,:] = data[nb_atoms*(Nrho+Nr)+interaction*Nr:nb_atoms*(Nrho+Nr)+interaction*Nr+Nr]
            interaction+=1
        return source,parameters, F,f,rep
    elif kind == "eam/fs":
        # reading 3 first comment lines as source for eam potential data
        source = eam[0].strip()+eam[1].strip()+eam[2].strip()
        # -- Parameters -- #
        atoms = eam[3].strip().split()[1:]
        nb_atoms = len(atoms)
        Nrho = int(eam[4].split()[0])       # Nrho (number of values for the embedding function F(rho))
        Nr = int(eam[4].split()[2])         # Nr (number of values for the effective charge function Z(r) and density function rho(r))
        drho = float(eam[4].split()[1])     # spacing in density space
        dr = float(eam[4].split()[3]) # spacing in distance space
        cutoff = float(eam[4].split()[4])
        atnumber,atmass,crystallatt,crystal = np.empty(nb_atoms,dtype=int),np.empty(nb_atoms),np.empty(nb_atoms),np.empty(nb_atoms).astype(np.str)
        for i in range(nb_atoms):
            # Fixme: The following lines assume that data occurs in blocks of
            # homogeneous width. This can break.
            l = len(eam[6].strip().split())
            row = int(5+i*((Nr*nb_atoms+Nrho)/l+1))
            atnumber[i] = int(eam[row].split()[0])
            atmass[i] = float(eam[row].split()[1])
            crystallatt[i] = float(eam[row].split()[2])
            crystal[i] = str(eam[row].split()[3])
        parameters = EAMParameters(atoms,atnumber,atmass,crystallatt,crystal,Nrho,Nr,drho,dr,cutoff)
        # -- Tabulated data -- #
        F,f,rep,data = np.empty((nb_atoms,Nrho)),np.empty((nb_atoms,nb_atoms,Nr)),np.empty((nb_atoms,nb_atoms,Nr)),np.empty(())
        eam = open(eam_file,'r')
        [eam.readline() for i in range(5)]
        for i in range(nb_atoms):
            eam.readline()
            data = np.append(data,np.fromfile(eam,count=Nrho+Nr*nb_atoms, sep=' '))
        data = np.append(data,np.fromfile(eam,count=-1, sep=' '))
        data = data[1:]
        for i in range(nb_atoms):
            F[i,:] = data[i*(Nrho+Nr*nb_atoms):Nrho+i*(Nrho+Nr*nb_atoms)]
            #f[i,:] = data[Nrho+i*(Nrho+Nr):Nrho+Nr+i*(Nrho+Nr)]
        interaction = 0
        for i in range(nb_atoms):
            for j in range(nb_atoms):
                f[i,j,:] = data[Nrho+j*Nr+i*(Nrho+Nr*nb_atoms):Nr+Nrho+j*Nr+i*(Nrho+Nr*nb_atoms)]
                if j < i :
                    rep[i,j,:] = data[nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr:nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr+Nr]
                    rep[j,i,:] = data[nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr:nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr+Nr]
                    interaction+=1
            rep[i,i,:] = data[nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr:nb_atoms*(Nrho+Nr*nb_atoms)+interaction*Nr+Nr]
            interaction+=1
        return source,parameters, F,f,rep
    else:
        print('Non supported eam file type')
        raise ValueError

