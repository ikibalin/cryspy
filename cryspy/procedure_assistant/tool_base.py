"""Exception and helpers shared by the assistant tool modules."""
from typing import List

import numpy

MAX_PARAMETERS_IN_LIST = 200


class ArgumentsNeeded(Exception):
    """Raised by a tool for the values to ask the user.

    l_missing: dicts with "name" (argument of the tool), "question" (what
    is asked, e.g. "the unit cell parameters") and optional JSON schema
    keys of the value ("type", "items", "minItems", "maxItems", "enum"),
    "hint" (format of the answer), "proposal" (a value known already, the
    user confirms it), "s_proposal" (its text)."""

    def __init__(self, l_missing: List[dict], message: str = ""):
        super().__init__(message)
        self.l_missing = l_missing
        self.message = message


def format_value(value) -> str:
    if isinstance(value, (float, numpy.floating)):
        return f"{float(value):.6g}"
    return str(value)
