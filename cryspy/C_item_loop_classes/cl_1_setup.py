from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Setup(ItemN):
    """
    Mandatory attributes:
        - wavelength

    Optional attributes:
        - field
        - offset_ttheta
        - offset_phi

    """
    ATTR_MANDATORY_NAMES = ("wavelength", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("wavelength", )

    ATTR_OPTIONAL_NAMES = ("field", "offset_ttheta", "offset_phi")
    ATTR_OPTIONAL_TYPES = (float, float, float)
    ATTR_OPTIONAL_CIF = ("field", "offset_2theta", "offset_phi")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("wavelength", "offset_ttheta", "offset_phi")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"offset_2theta":0.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "setup"

    def __init__(self, **kwargs) -> NoReturn:
        super(Setup, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"wavelength": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class SetupL(LoopN):
    """
    Description of Setup in loop.

    """
    ITEM_CLASS = Setup
    ATTR_INDEX = None
    def __init__(self, loop_name = None) -> NoReturn:
        super(SetupL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#   loop_
#  _setup_wavelength   0.84
#  _setup_field        1.00
#  _setup_offset_2theta -0.385
#  _setup_offset_phi -0.385  
#     0.84 1.0 -0.385 -0.385
#     1.0  1.5 0.7 0.3
#   """

# obj = SetupL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj[0], end="\n\n")
