import numpy as np

def random_three_vector():
    """Generates a random 3D unit vector (direction) with a
    uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    """

    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)


def random_velocities(individual):
    velocities = np.zeros((len(individual), 3), dtype=np.float)
    for row in velocities:
        vx, vy, vz = random_three_vector()
        row[0], row[1], row[2] = vx, vy, vz
    individual.set_velocities(velocities)

