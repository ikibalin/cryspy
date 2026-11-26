"""Exceptions and warnings raised by CrysPy.

Catch every CrysPy error with ``except cryspy.CrysPyError``. Each error
also derives from the built-in exception that was raised at the same
place before, so existing ``except ValueError`` / ``except AttributeError``
/ ``except FileNotFoundError`` code keeps working.
"""


from typing import List


class CrysPyError(Exception):
    """Base class of all CrysPy errors."""


class ArgumentsNeeded(Exception):
    """Raised by a function for the values to ask the user.

    A procedure, a method or a tool of AI assistants raises it when it is
    in doubt about a value (several possible values, a guess the user
    should confirm); a program shows the questions and calls the function
    again with the answers as arguments of the same names.

    l_missing: dicts with "name" (argument of the function), "question"
    (what is asked, e.g. "the unit cell parameters") and optional JSON
    schema keys of the value ("type", "items", "minItems", "maxItems",
    "enum"), "hint" (format of the answer), "options" (possible answers,
    shown as buttons; the user may type another one), "proposal" (a value
    known already, the user confirms it), "s_proposal" (its text)."""

    def __init__(self, l_missing: List[dict], message: str = ""):
        super().__init__(message)
        self.l_missing = l_missing
        self.message = message


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
