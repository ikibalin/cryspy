"""Classes Exclude, ExcludeL."""
from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Exclude(ItemN):
    """
    Exclued class.

    Attributes
    ----------
        - ttheta_min (mandatory)
        - ttheta_max (mandatory)
        - id
        - phi_min
        - phi_max
    """

    ATTR_MANDATORY_NAMES = ("ttheta_min", "ttheta_max")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("2theta_min", "2theta_max")

    ATTR_OPTIONAL_NAMES = ("id", "phi_min", "phi_max")
    ATTR_OPTIONAL_TYPES = (str, float, float)
    ATTR_OPTIONAL_CIF = ("id", "phi_min", "phi_max")

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

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "exclude"

    def __init__(self, **kwargs) -> NoReturn:
        super(Exclude, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"ttheta_min": 0., "ttheta_max": 0., "phi_min": -90.,
                 "phi_max": -90.}

        # defined for ani integer and float parameters
        D_MAX = {"ttheta_min": 180., "ttheta_max": 180., "phi_min": 90.,
                 "phi_max": 90.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class ExcludeL(LoopN):
    """Description of Exclude in loop."""

    ITEM_CLASS = Exclude
    ATTR_INDEX = "id"

    def __init__(self, loop_name=None) -> NoReturn:
        super(ExcludeL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#  loop_
#  _exclude_id
#  _exclude_ttheta_min
#  _exclude_ttheta_max
#   1   4.0  12.0
#   2  30.0  45.0
#   3  58.0  63.0
# """

# obj = ExcludeL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["2"], end="\n\n")
# print(obj.ttheta_min, end="\n\n")
