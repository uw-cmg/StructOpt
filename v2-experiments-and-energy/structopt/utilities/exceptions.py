class UnknownState(Exception):
    pass

class Running(Exception):
    pass

class Queued(Exception):
    pass

class Submitted(Exception):
    def __init__(self, jobdir):
        self.jobdir = jobdir
    def __str__(self):
        return repr(self.jobdir)
