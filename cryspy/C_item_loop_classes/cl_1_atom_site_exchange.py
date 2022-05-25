"""Describe AtomSiteSusceptibility, AtomSiteSusceptibilityL."""
from typing import NoReturn
import math
import numpy

from cryspy.A_functions_base.function_1_algebra import calc_m_sigma
from cryspy.A_functions_base.function_1_atomic_vibrations import \
    vibration_constraints
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

import cryspy.A_functions_base.local_susceptibility as local_susceptibility

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

na = numpy.newaxis


class AtomSiteExchange(ItemN):
    """Exchange interaction (only in case of G-tensor)
    """

    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = (
        "j_type", "j_11", "j_22", "j_33", "j_12",
        "j_13", "j_23", )
    ATTR_OPTIONAL_TYPES = (str, float, float, float, float, float, float, )
    ATTR_OPTIONAL_CIF = (
        "J_type", "J_11", "J_22", "J_33", "J_12",
        "J_13", "J_23", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("j_11", "j_22", "j_33", "j_12", "j_13", "j_23",)
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"j_11": "{:.5f}", "j_22": "{:.5f}", "j_33": "{:.5f}",
                 "j_12": "{:.5f}", "j_13": "{:.5f}", "j_23": "{:.5f}",}

    # constraints on the parameters
    D_CONSTRAINTS = {"j_type": ["Jiso", ],
                     }

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_exchange"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteExchange, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


    def apply_j_iso_constraint(self):
        """Isotropic constraint on exchange."""
        if not(self.is_attribute("j_type")):
            return
        j_type = self.j_type
        if j_type.lower().startswith("jiso"):
            self.__dict__["j_22"] = self.j_11
            self.__dict__["j_33"] = self.j_11
            self.__dict__["j_12"] = 0. 
            self.__dict__["j_13"] = 0. 
            self.__dict__["j_23"] = 0. 
            self.__dict__["j_22_sigma"] = self.j_11_sigma
            self.__dict__["j_33_sigma"] = self.j_11_sigma
            self.__dict__["j_12_sigma"] = 0. 
            self.__dict__["j_13_sigma"] = 0. 
            self.__dict__["j_23_sigma"] = 0. 
            self.__dict__["j_22_refinement"] = False
            self.__dict__["j_33_refinement"] = False
            self.__dict__["j_12_refinement"] = False
            self.__dict__["j_13_refinement"] = False
            self.__dict__["j_23_refinement"] = False
            self.__dict__["j_22_constraint"] = True
            self.__dict__["j_33_constraint"] = True
            self.__dict__["j_12_constraint"] = True
            self.__dict__["j_13_constraint"] = True
            self.__dict__["j_23_constraint"] = True

    def get_flags_exchange(self):
        res = numpy.array([
            self.j_11_refinement, self.j_22_refinement, self.j_33_refinement,
            self.j_12_refinement, self.j_13_refinement, self.j_23_refinement], dtype=bool)
        return res



class AtomSiteExchangeL(LoopN):
    """Exchange parameters for magnetic atoms.

    Methods
    -------
        - apply_space_group_constraint
        - apply_chi_iso_constraint
    """

    ITEM_CLASS = AtomSiteExchange
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteExchangeL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def apply_j_iso_constraint(self):
        """Apply isotropic constraint on susceptibility."""
        for item in self.items:
            item.apply_j_iso_constraint()

    def get_flags_exchange(self):
        flags_exchange = numpy.stack([item.get_flags_exchange() for item in self.items], axis=0).transpose()
        return flags_exchange

