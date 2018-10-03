import functools
import weakref


def inner_lazy(method, is_property):
    values = weakref.WeakKeyDictionary()

    @functools.wraps(method)
    def wrapper(self):
        try:
            return values[self]
        except KeyError:
            values[self] = method(self)
            return values[self]

    def _lazy_reset(self):
        try:
            del values[self]
        except KeyError:
            pass
    wrapper._lazy_reset = _lazy_reset

    if is_property:
        wrapper = property(wrapper)

    return wrapper


def lazyproperty(method):
    return inner_lazy(method, True)


def lazy(is_property):
    if callable(is_property):
        # Then is_property was not defined and we should use the default value
        # Moreover, is_property is actually the method
        return inner_lazy(is_property, False)
    else:
        # is_property is a bool
        def wrapper(method):
            return inner_lazy(method, is_property)
        return wrapper

