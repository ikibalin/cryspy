"""DensityPoint and DensityPointL classes."""
from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

import numpy

class ChannelPlusMinus(ItemN):
    """
    Magnetization density in channel plus minus.

    Attributes
    ----------
        - numerator_x, numerator_y, numerator_z, denominator_xyz (mandatory)
        - index_x, index_y, index_z (optional)
        - density_plus
        - density_minus
        - point_multiplicity (optional)
    """

    ATTR_MANDATORY_NAMES = ("numerator_x", "numerator_y", "numerator_z", "denominator_xyz",
        "density_plus", "density_minus", )
    ATTR_MANDATORY_TYPES = (int, int, int, int,
        float, float, )
    ATTR_MANDATORY_CIF = ("numerator_x", "numerator_y", "numerator_z", "denominator_xyz",
        "density_plus", "density_minus", )

    ATTR_OPTIONAL_NAMES = ("index_x", "index_y", "index_z", "point_multiplicity", )
    ATTR_OPTIONAL_TYPES = (int, int, int, int)
    ATTR_OPTIONAL_CIF = ("index_x", "index_y", "index_z", "point_multiplicity", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"density_plus": "{:.5f}", "density_minus": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "channel_plus_minus"

    def __init__(self, **kwargs) -> NoReturn:
        super(ChannelPlusMinus, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"density_plus": 0., "density_minus": 0.,
                 "point_multiplicity": 0,}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class ChannelPlusMinusL(LoopN):
    """
    Magnetization density in channel chi in loop

    Attributes
    ----------
        - numerator_x, numerator_y, numerator_z, denominator_xyz (mandatory)
        - index_x, index_y, index_z (optional)
        - density
        - chi_11, chi_12, chi_13, chi_21, chi_22, chi_23, chi_31, chi_32, chi_33
        - atom_multiplicity
        - point_multiplicity (optional)
    """

    ITEM_CLASS = ChannelPlusMinus
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ChannelPlusMinusL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

