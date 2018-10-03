import sys
import functools


def get_rank():
    if 'mpi4py' in sys.modules:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_rank()
    else:
        return 0


def get_size():
    if 'mpi4py' in sys.modules:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_size()
    else:
        return 0


def root(method=None, broadcast=True):
    """A decorator to make the function only run on the root node. The returned
    data from the root is then broadcast to all the other nodes and each node
    returns the root's data.
    """
    # This is a complicated decorator. See the following link for some context:
    # https://blogs.it.ox.ac.uk/inapickle/2012/01/05/python-decorators-with-optional-arguments/
    if method is None:
        return functools.partial(root, broadcast=broadcast)

    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        if broadcast and 'mpi4py' in sys.modules:
            from mpi4py import MPI
            if MPI.COMM_WORLD.Get_rank() == 0:
                data = method(*args, **kwargs)
            else:
                data = None
            if hasattr(data, 'bcast'):
                data.bcast()
            else:
                data = MPI.COMM_WORLD.bcast(data, root=0)
        else:
            data = method(*args, **kwargs)
        return data
    wrapper.__doc__ += "\n\n(@root) Designed to run on the root node only.\n"
    return wrapper


def parallel(method):
    """A decorator that does nothing except document that the function is
    designed to run in parallel.
    """
    @functools.wraps(method)
    def wrapper(*args, **kwargs):
        return method(*args, **kwargs)
    wrapper.__doc__ += ("\n\n(@parallel) Designed to run code that runs differently on different cores.\n"
                        "The MPI functionality should be implemented inside these functions.\n")
    return wrapper


def single_core(method):
    """A place holder decorator that does nothing except document that the function is designed to be run on a single core."""
    return method


def allgather(stuff, stuffs_per_core):
    """Performs an MPI.allgather on a selection of data and uses stuffs_per_core
    to parse out the correct information and return it.

    Args:
        stuff (any): any piece of data (e.g. fitnesses), some of which have been updated
            on their respective cores and some of which haven't. each piece of data should
            be of the same length
        stuffs_per_core (dict<int, list<int>>): a dictionary containing a mapping of the
            cores that contain the correct information to the corresponding indices in the
            pieces of data

    Returns:
        type(stuff): the correct stuff is returned on each core

    Example:
        This is going to take the values::

            values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        and convert each of them to strings.

        In this example there are 5 cores, so ``stuffs_per_core`` looks like::

            stuffs_per_core = {0: [0, 5], 1: [1, 6], 2: [2, 7], 3: [3, 8], 4: [4, 9]}

        Now for the code that precedes ``allgather()`` and then calls ``allgather()``::

            # This for-loop modifies different parts of `values` on each core by
            # converting some elements in `values` from an int to a str.
            # We then want to collect the values that each core independently updated
            # and allgather them so that every core has all of the updated values,
            # even though each core only did part of the work.
            for i in stuffs_per_core[rank]:
                values[i] = str(inds[i])
            x = allgather(values, stuffs_per_core)
            print(x)  # returns:  ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    """
    if hasattr(stuff, 'allgather'):
        raise TypeError('instance of {} has an `allgather` function that should be used instead'.format(stuff.__class__.__name__))

    from mpi4py import MPI
    # The lists in stuffs_per_core all need to be of the same length
    max_stuffs_per_core = max(len(stuffs) for stuffs in stuffs_per_core.values())
    for rank, stuffs in stuffs_per_core.items():
        while len(stuffs) < max_stuffs_per_core:
            stuffs.append(None)

    all_stuffs_per_core = MPI.COMM_WORLD.allgather(stuff)
    correct_stuff = [None for _ in range(len(stuff))]
    for rank, indices in stuffs_per_core.items():
        for index in indices:
            if index is not None:
                correct_stuff[index] = all_stuffs_per_core[rank][index]

    # If something didn't get sent, use the value on the core
    for i, thing in enumerate(correct_stuff):
        if thing is None:
            correct_stuff[i] = stuff[i]

    return correct_stuff


def parse_MPMD_cores_per_structure(value):
    """Converts an input ``value`` from a value in the parameter file into a ``{'min': ..., 'max': ...}`` dictionary."""
    if isinstance(value, int):
        if int == 0:
            return None
        else:
            return {'min': value, 'max': value}
    elif isinstance(value, str):
        if value == 'any':
            from mpi4py import MPI
            return {'min': 1, 'max': MPI.COMM_WORLD.Get_size()}
        elif '-' not in value:
            value = int(value)
            return {'min': value, 'max': value}
        else:
            min_, max_ = value.split('-')
            return {'min': int(min_), 'max': int(max_)}
    else:
        raise TypeError("'MPMD_cores_per_structure' must be an 'int' or 'str'")


class SingleCorePerIndividual(object):
    """
    This is an untested context manager for using the mpi4py parallelization technique.
    I'm not sure we should use it, but I'm thinking on it.
    The nice thing is that it is self-contained reasonably well.
    The bad thing is that the user has to manually set fields on the individuals before
    the context manager exits in order for it to work properly,
    and there isn't a great way to error check whether they did that or not.

    Example usage:

        def calculate_energies(population):
            to_fit = [individual for individual in population if individual.LAMMPS is None]
            with SingleCorePerIndividual(population, to_run, "LAMMPS") as to_compute:
                for individual in to_compute:
                    individual.LAMMPS = individual.fitnesses.LAMMPS.fitness(individual)
            return [individual.LAMMPS for individual in population]
    """
    def __init__(self, population, to_run, attr_name):
        import gparameters
        self.population = population
        self.to_run = to_run
        self.rank = gparameters.mpi.rank
        self.ncores = gparameters.mpi.ncores

    def __enter__(self):
        self.individuals_per_core = {r: [] for r in range(self.ncores)}
        for i, individual in enumerate(self.to_run):
            self.individuals_per_core[i % self.ncores].append(individual)
        return self.individuals_per_core[self.rank]

    def __exit__(self, exception_type, exception_value, traceback):
        positions_per_core = {rank: [self.population.position(individual) for individual in individuals] for rank, individuals in self.individuals_per_core.items()}
        self.allgather([getattr(individual, self.attr_name) for individual in self.population])
        allgather(self.results, positions_per_core)


"""
@root
def fitness(population, parameters):
    from mpi4py import MPI

    to_fit = [individual for individual in population if not individual._fitted]

    if to_fit:
        spawn_args = [individual.fitnesses.FEMSIM.get_spawn_args(individual) for individual in to_fit]
        MPMD(to_fit, spawn_args, parameters)

        # Collect the results for each chisq and return them
        logger = logging.getLogger('output')
        for i, individual in enumerate(to_fit):
            vk = individual.fitnesses.FEMSIM.get_vk_data()
            individual.FEMSIM = individual.fitnesses.FEMSIM.chi2(vk)
            logger.info('Individual {0} for FEMSIM evaluation had chisq {1}'.format(i, individual.FEMSIM))

    return [individual.FEMSIM for individual in population]
"""

@root(broadcast=False)
def MPMD(self, to_run, spawn_args, parameters):
    """This is a function to run MPMD, with an example for FEMSIM above.
    I'm not confident we should use it yet, but only because I haven't tested it; I'm relatively confident it will work.
    """
    ncores = parameters.ncores
    cores_per_individual = ncores // len(to_run)
    # Round cores_per_individual down to nearest power of 2
    if cores_per_individual == 0:
        # We have more individuals than cores, so each fitness scheme needs to be run multiple times
        cores_per_individual = 1

    pow(2.0, math.floor(math.log2(cores_per_individual)))
    minmax = parse_MPMD_cores_per_structure(parameters.MPMD)
    assert minmax['max'] >= minmax['min'] > 0
    if cores_per_individual < minmax['min']:
        cores_per_individual = minmax['min']
    elif cores_per_individual > minmax['max']:
        cores_per_individual = minmax['max']

    # Setup each individual and get the inputs for each individual that need to be passed into the spawn
    multiple_spawn_args = defaultdict(list)
    for individual in to_run:
        for arg_name, arg in spawn_args.items():
            multiple_spawn_args[arg_name].append(arg)

    # Make sure each key in `multiple_spawn_args` has the same number of elements and set `count` to that value
    count = -1
    for key, value in multiple_spawn_args.items():
        if count != -1:
            assert len(value) == count
        count = len(value)

    # Create MPI.Info objects from the kwargs dicts in multiple_spawn_args['info'] for each rank
    infos = [MPI.Info.Create() for _ in multiple_spawn_args['info']]
    for i, info in enumerate(infos):
        for key, value in multiple_spawn_args['info'][i].items():
            info.Set(key, value)

    # Run the multiple spawn
    individuals_per_iteration = ncores // cores_per_individual
    individuals_per_iteration = min(individuals_per_iteration, len(to_run))
    num_iterations = math.ceil(len(to_run) / individuals_per_iteration)
    for i in range(num_iterations):
        j = i * individuals_per_iteration
        individuals_this_iteration = individuals_per_iteration
        if i == len(to_run) // individuals_per_iteration:  # The last iteration may not be exactly individuals_per_iteration
            individuals_this_iteration = len(to_run) % individuals_per_iteration
            cores_per_individual = ncores // individuals_this_iteration # All the cores available should be used
        print("Spawning {} femsim processes, each with {} cores".format(individuals_this_iteration, cores_per_individual))
        intercomm = MPI.COMM_SELF.Spawn_multiple(command=multiple_spawn_args['command'][j:j+individuals_this_iteration],
                                                 args=multiple_spawn_args['args'][j:j+individuals_this_iteration],
                                                 maxprocs=[cores_per_individual]*individuals_this_iteration,
                                                 info=infos[j:j+individuals_this_iteration]
                                                 )
        # Disconnect the child processes
        intercomm.Disconnect()

