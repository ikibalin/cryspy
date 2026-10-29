"""Marking of public functions and methods planned for removal."""
import functools
import warnings


def deprecated(func):
    """Warn with DeprecationWarning on each call of the decorated function.

    Used for public methods that nothing in CrysPy, CrysPy Editor and
    PyColl5 calls any more: user scripts get one release to move away
    before the method is removed.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            f"{func.__qualname__} is deprecated and will be removed in "
            "the next CrysPy release.", DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    wrapper.__deprecated__ = True  # left out of the methods for AI assistants
    return wrapper
