import json
import random


class Random(object):
    def __init__(self):
        self.step_number = 0

        self.args = self.read_inputs()
        random.seed(self.args.get('seed', None))


    def read_inputs(self):
        args = json.load(open('inputs.json'))
        return args


    def update_parameters(self, **kwargs):
        self.step_number += 1
        for key, value in kwargs.items():
            self.args[key] = value


    def fitness(self, individuals):
        return [random.random() for _ in individuals]
