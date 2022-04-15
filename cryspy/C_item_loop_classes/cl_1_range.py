from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Range(ItemN):
    """Contains range information.
    
    Attributes
    ----------
        - 
    """
    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("ttheta_min", "ttheta_max", "phi_min", "phi_max",
                           "time_min", "time_max", "gamma_min", "gamma_max", "nu_min", "nu_max")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("2theta_min", "2theta_max", "phi_min", "phi_max",
                         "time_min", "time_max", "gamma_min", "gamma_max", "nu_min", "nu_max")

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

    PREFIX = "range"

    def __init__(self, **kwargs) -> NoReturn:
        super(Range, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"ttheta_min": 0., "ttheta_max": 0., "phi_min": -90.,
                 "phi_max": -90., "time_min": 0., "time_max": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"ttheta_min": 180., "ttheta_max": 180., "phi_min": 90.,
                 "phi_max": 90.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)



class RangeL(LoopN):
    """
    Description of Range in loop.

    """
    ITEM_CLASS = Range
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(RangeL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#   loop_
#  _range_phi_min
#  _range_2theta_min
#  _range_2theta_max
#  7 4 80
#  7 -9 780
#   """

# obj = RangeL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj[0], end="\n\n")
