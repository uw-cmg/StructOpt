import sys
import random

import numpy as np

import structopt
from structopt.common.population import Population

parameters = structopt.setup('structopt.in.json')
random.seed(parameters.seed)
np.random.seed(parameters.seed)

population = Population(parameters=parameters)
