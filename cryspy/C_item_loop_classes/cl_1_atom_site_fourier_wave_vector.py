"""AtomSiteFourierWaveVector, AtomSiteFourierWaveVectorL classes are given (mCIF category atom_site_Fourier_wave_vector)."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomSiteFourierWaveVector(ItemN):
    """Wave vector of a Fourier term of the modulations (msCIF).

``q1_coeff`` is the harmonic of the propagation vector (only 1 is
supported); ``x``, ``y``, ``z`` give the vector itself.
    """

    ATTR_MANDATORY_NAMES = ("seq_id", )
    ATTR_MANDATORY_TYPES = (int, )
    ATTR_MANDATORY_CIF = ("seq_id", )

    ATTR_OPTIONAL_NAMES = ("q1_coeff", "x", "y", "z", )
    ATTR_OPTIONAL_TYPES = (int, float, float, float, )
    ATTR_OPTIONAL_CIF = ("q1_coeff", "x", "y", "z", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {}

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_fourier_wave_vector"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteFourierWaveVector, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class AtomSiteFourierWaveVectorL(LoopN):
    """Wave vectors of the Fourier terms of the modulations (msCIF).
    """

    ITEM_CLASS = AtomSiteFourierWaveVector
    ATTR_INDEX = "seq_id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteFourierWaveVectorL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
