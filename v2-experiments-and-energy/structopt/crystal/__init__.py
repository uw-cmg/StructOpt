from structopt.common.individual import Individual


class Crystal(Individual):
    """A stucture with periodic boundary conditions."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.set_pbc([True, True, True])

