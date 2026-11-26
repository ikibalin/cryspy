"""Exception and helpers shared by the assistant tool modules."""

import numpy

from cryspy.exceptions import ArgumentsNeeded  # noqa: F401 (imported from here by the tools)

MAX_PARAMETERS_IN_LIST = 200


def format_value(value) -> str:
    if isinstance(value, (float, numpy.floating)):
        return f"{float(value):.6g}"
    return str(value)
