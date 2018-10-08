import sys
import random
import logging
import numpy as np
from scipy.special import sph_harm
from scipy.linalg import norm
import structopt
from structopt.common.population import Population

def update_particle(individual, best_swarm, best_particle, omega, phi_p, phi_g):
    natoms = len(individual.positions)
    """
    dist = 0.0
    max_dist = 0.0
    max_individual = None
    thresh = 1e-3
    old_individual = individual.copy()
    itr = 0
    while dist < thresh and itr < 1:
        rp = np.random.rand(natoms,3)
        rg = np.random.rand(natoms,3)
        
        individual.set_velocities( omega * old_individual.velocities +
                phi_p * rp * (best_particle.positions - old_individual.positions) +
                phi_g * rg * (best_swarm.positions - old_individual.positions) )
        individual.set_positions(old_individual.positions + individual.velocities)
        individual._Q_l = np.array([])
        dist = distance_BCM(old_individual, individual, 3.0)
        if max_dist < dist:
            max_dist = dist
            max_individual = individual.copy()
        itr += 1
    individual = max_individual.copy()
    #print("Individual %s: %s"%(individual.id, max_dist))
    """
    rp = np.random.rand(natoms,3)
    rg = np.random.rand(natoms,3)
    #choice = np.array([[random.choice([1, 0])] for i in range(natoms)])
    individual.set_velocities( omega * individual.velocities +
            phi_p * rp * (best_particle.positions - individual.positions) +
            phi_g * rg * (best_swarm.positions - individual.positions) )
    individual.set_positions(individual.positions + individual.velocities)
    return None

def distance_BCM(individualA, individualB, cutoff=3.0):
    dist = 0
    l_set = range(2, 14, 2)
    set_Q_l(individualA, l_set, cutoff)
    set_Q_l(individualB, l_set, cutoff)
    dist = sum( (individualA._Q_l - individualB._Q_l)**2 )
    return dist/len(l_set)

def set_Q_l(individual, l_set, cutoff=3.0):
    if len(individual._Q_l) > 0:
        return
    xyz_list, r_list = get_bonds(individual, cutoff)
    tan_theta = xyz_list[:,1] / xyz_list[:,0] # y/x
    tan_theta[np.isnan(tan_theta)] = 0
    theta = np.arctan(tan_theta)
    phi = np.arccos(xyz_list[:,2]/r_list) # z/r
    myQ_l = np.array([
             ((4*np.pi/(2*l+1)) * 
                    sum([abs(ave_Q_ml(m, l, theta, phi))**2 
                for m in range(-l, l+1)])) ** 0.5 
             for l in l_set])
    individual._Q_l = myQ_l
    return 

def get_bonds(individual, cutoff=3.0):
    natoms = len(individual.positions)
    #xyz_list = np.array([ individual.positions[i] - individual.positions[j] 
    #            for i in range(natoms-1) for j in range(i+1, natoms)] )
    xyz_list = np.repeat(individual.positions, range(natoms-1,-1,-1), axis=0) - np.vstack([individual.positions[j:,:] for j in range(1,natoms)])
    r_list = norm(xyz_list, axis=1)
    xyz_list = xyz_list[r_list < cutoff]
    r_list = r_list[r_list < cutoff]
    return xyz_list, r_list

def ave_Q_ml(m, l, theta, phi):
    cnt = len(theta)
    return np.sum(sph_harm(m, l, theta, phi))/cnt if cnt else 0


