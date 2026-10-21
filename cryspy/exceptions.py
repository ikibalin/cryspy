"""Exceptions and warnings raised by CrysPy.

Catch every CrysPy error with ``except cryspy.CrysPyError``. Each error
also derives from the built-in exception that was raised at the same
place before, so existing ``except ValueError`` / ``except AttributeError``
/ ``except FileNotFoundError`` code keeps working.
"""


class CrysPyError(Exception):
    """Base class of all CrysPy errors."""


class RcifFormatError(CrysPyError, ValueError):
    """The rcif text, or a value in it, cannot be understood.

    The message names the tag and the value when they are known.
    """


class RcifFileNotFoundError(CrysPyError, FileNotFoundError):
    """The rcif file to read does not exist."""


class RefinementError(CrysPyError, AttributeError):
    """A procedure (RhoChi refinement, MEM reconstruction) cannot run on
    the given object, e.g. there is no crystal or no experiment in it."""


class CrysPyWarning(UserWarning):
    """Something in the input is ignored or replaced; the calculation goes on."""
