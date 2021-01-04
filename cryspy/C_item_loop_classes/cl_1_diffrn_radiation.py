from typing import NoReturn

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class DiffrnRadiation(ItemN):
    """
    Polarization of the neutron beam.

    Attributes:
        - polarization  (float)
        - efficiency (float)

    """
    ATTR_MANDATORY_NAMES = ("polarization", "efficiency")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("polarization", "efficiency")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("polarization", "efficiency")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"polarization": 1., "efficiency": 1.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "diffrn_radiation"

    def __init__(self, **kwargs) -> NoReturn:
        super(DiffrnRadiation, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"polarization": -1., "efficiency": -1.}

        # defined for ani integer and float parameters
        D_MAX = {"polarization": 1., "efficiency": 1.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class DiffrnRadiationL(LoopN):
    """
    Description of beam polarization in loop.

    """
    ITEM_CLASS = DiffrnRadiation
    ATTR_INDEX = None
    def __init__(self, loop_name = None) -> NoReturn:
        super(DiffrnRadiationL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
# loop_
#  _diffrn_radiation_polarization 
#  _diffrn_radiation_efficiency   
#  1. 1.
#  1. -0.87()
#  5() -7
# """

# obj = DiffrnRadiationL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.get_variable_names(), end="\n\n")
