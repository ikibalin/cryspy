"""CellWaveVector, CellWaveVectorL classes are given (mCIF category cell_wave_vector)."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class CellWaveVector(ItemN):
    """Propagation (modulation) vector k of the magnetic structure.

The components x, y, z are in reciprocal lattice units of the cell
(mCIF/msCIF ``_cell_wave_vector_x``). Only one vector is supported
(3+1 dimensional superspace groups). Standard uncertainties such as
``0.0197(2)`` are read; the refinement of k is not implemented.
    """

    ATTR_MANDATORY_NAMES = ("seq_id", "x", "y", "z", )
    ATTR_MANDATORY_TYPES = (int, float, float, float, )
    ATTR_MANDATORY_CIF = ("seq_id", "x", "y", "z", )

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("x", "y", "z")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"x": "{:.5f}", "y": "{:.5f}", "z": "{:.5f}"}

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "cell_wave_vector"

    def __init__(self, **kwargs) -> NoReturn:
        super(CellWaveVector, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class CellWaveVectorL(LoopN):
    """Propagation vectors of the magnetic structure (one is supported).
    """

    ITEM_CLASS = CellWaveVector
    ATTR_INDEX = "seq_id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(CellWaveVectorL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
