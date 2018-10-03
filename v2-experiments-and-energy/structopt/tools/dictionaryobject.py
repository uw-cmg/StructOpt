import os.path
from collections import Mapping


class DictionaryObject(dict):
    """A dictionary-like object that allows attribute access (both getting and setting) via an object's dot notation."""
    # Unless you really like python, you probably don't want to worry about this implementation.
    def __init__(self, dictionary):
        for k, v in dictionary.items():
            setattr(self, k, v)

    @staticmethod
    def _render(obj):
        # JSON can only handle a few input types: int/float, str, bool, null, dict, and list
        if isinstance(obj, int) or isinstance(obj, float) or isinstance(obj, bool) or obj is None:
            return obj
        if isinstance(obj, str):
            return os.path.expandvars(obj)
        elif isinstance(obj, dict):
            return DictionaryObject(obj)
        elif isinstance(obj, list):
            return [DictionaryObject._render(v) for v in obj]
        raise TypeError("Unknown type given to DictionaryObject._render: {}".format(type(obj)))

    def __getattr__(self, key):
        return super().__getitem__(key)

    __detattr__ = dict.__delitem__

    def __setattr__(self, key, value):
        super().__setitem__(key, DictionaryObject._render(value))

    def __setitem__(self, key, value):
        super().__setitem__(key, DictionaryObject._render(value))

    def setdefault(self, key, default):
        if key not in self:
            self[key] = default

    def update(self, other=None, **kwargs):
        if other is not None:
            for k, v in other.items() if isinstance(other, Mapping) else other:
                self[k] = v
        for k, v in kwargs.items():
            self[k] = v

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        # Restore instance attributes
        self.__dict__.update(state)

