import os
try:
    from ase import Atom, Atoms
    from ase.optimize import BFGS
    from ase.units import GPa
    from ase.calculators.neighborlist import NeighborList
except ImportError:
    pass
from StructOpt.inp_out.write_xyz import write_xyz
from StructOpt.tools.find_defects import find_defects
from StructOpt.tools.check_cell_type import check_cell_type
from StructOpt.fingerprinting import get_fingerprint
from StructOpt.tools.lammps import LAMMPS
from StructOpt.tools.setup_energy_calculator import setup_energy_calculator
import numpy
import math
import json
try:
    from mpi4py import MPI
except ImportError:
    pass
import logging
import pdb
import shutil
import time
from scipy import spatial
import random

def check_structures(Optimizer,indiv):
        # Check for atoms that are too close or out of constrained location
        stro = ''
        for ind in range(len(indiv)):
            totalsol = indiv[ind][0].copy()
            if Optimizer.constrain_position:
                if Optimizer.structure=='Defect':
                    totalsol, stro = constrain_positions(totalsol, Optimizer.solidbulk, Optimizer.sf)
            min_len=0.7
            if not Optimizer.fixed_region:
                totalsol, stro = check_min_dist(Optimizer, totalsol, Optimizer.structure, totalsol.get_number_of_atoms(), min_len, stro) 
            indiv[ind][0] = totalsol.copy()
        return indiv, stro

def constrain_positions(indiv, bulk, sf):
        STR=''
        ts = indiv.copy()
        indc,indb,vacant,swap,stro = find_defects(ts,bulk,0)
        sbulk = bulk.copy()
        bcom = sbulk.get_center_of_mass()
        com = indc.get_center_of_mass()
        dist = (sum((bcom[i] - com[i])**2 for i in range(3)))**0.5
        if dist > sf:
            STR+='Shifting structure to within region\n'
            r = random.random()*sf
            comv = numpy.linalg.norm(com)
            ncom = [one*r/comv for one in com]
            trans = [ncom[i]-com[i] for i in range(3)]
            indices = []
            for one in indc:
                id = [atm.index for atm in totalsol if atm.position[0]==one.position[0] and atm.position[1]==one.position[1] and atm.position[2]==one.position[2]][0]
                ts[id].position += trans
        return ts, STR
    
def check_min_dist(Optimizer, totalsol, type='Defect', nat=None, min_len=0.7, STR=''):
        if type=='Defect' or type=='Crystal' or type=='Surface':
            if nat==None:
                nat=len(totalsol)
            cutoffs=[2.0 for one in totalsol]
            nl=NeighborList(cutoffs,bothways=True,self_interaction=False)
            nl.update(totalsol)
            for one in totalsol[0:nat]:
                nbatoms=Atoms()
                nbatoms.append(one)
                indices, offsets=nl.get_neighbors(one.index)
                for index, d in zip(indices,offsets):
                    index = int(index)
                    sym=totalsol[index].symbol
                    pos=totalsol[index].position + numpy.dot(d,totalsol.get_cell())
                    at=Atom(symbol=sym,position=pos)
                    nbatoms.append(at)
                while True:
                    dflag=False
                    for i in range(1,len(nbatoms)):
                        d=nbatoms.get_distance(0,i)
                        if d < min_len:
                            nbatoms.set_distance(0,i,min_len+.01,fix=0.5)
                            STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
                            dflag=True
                    if dflag==False:
                        break
                for i in range(len(indices)):
                    totalsol[indices[i]].position=nbatoms[i+1].position
                totalsol[one.index].position=nbatoms[0].position
                nl.update(totalsol)
        elif type=='Cluster':
            if not 'LAMMPS' in Optimizer.modules:
                for i in range(len(totalsol)):
                    for j in range(len(totalsol)):
                        if i != j:
                            d=totalsol.get_distance(i,j)
                            if d < min_len:
                                totalsol.set_distance(i,j,min_len,fix=0.5)
                                STR+='--- WARNING: Atoms too close (<0.7A) - Implement Move ---\n'
            else:
                rank = MPI.COMM_WORLD.Get_rank()
                #logger = logging.getLogger(Optimizer.loggername)
                R = totalsol.arrays['positions']
                tol = 0.01
                epsilon = 0.05
                fix = 0.5
                
                closelist = numpy.arange(len(totalsol))
                iter = 0
                while len(closelist) > 0 and iter<2:
                    iter+=1
                    closelist = []
                    dist=spatial.distance.cdist(R,R)
                    numpy.fill_diagonal(dist,1.0)
                    smalldist = numpy.where(dist < min_len-tol)
   
                    for ind in range(len(smalldist[0])):
                        i = smalldist[0][ind]
                        j = smalldist[1][ind]
                        if i < j and dist[i][j] < min_len-tol:
                            closelist.append(i)
                            closelist.append(j)
                            if dist[i][j] > epsilon:
                                x = 1.0 - min_len / dist[i][j]
                                D = R[j]-R[i]
                                R[i] += (x * fix) * D
                                R[j] -= (x * (1.0 - fix)) * D
                            else:
                                R[i] += [0.2, 0.0, 0.0]
                                R[j] -= [0.2, 0.0, 0.0]
                            R2P = [R[i],R[j]]
                            dist2P=spatial.distance.cdist(R2P,R)
                            dist[i] = dist2P[0]
                            dist[j] = dist2P[1]
                            for k in range(len(R)):
                                dist[k][i] = dist[i][k]
                                dist[k][j] = dist[j][k]
                    closelist=list(set(closelist))
                    closelist.sort()
                    #if len(closelist) != 0:
                    #    logger.info('M:iter {0}, closelist size {1}'.format(iter,len(closelist)))
        else:
            print 'WARNING: In Check_Min_Dist in EvalEnergy: Structure Type not recognized'
        return totalsol, STR
