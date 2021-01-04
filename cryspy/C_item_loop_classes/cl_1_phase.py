"""Description of Phase and PhaseL classes."""
from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Phase(ItemN):
    """Phase class.

    Attributes
    ----------
        - label (mandatory)
        - scale, igsize, u, v, w, x, y (optional)
    """

    ATTR_MANDATORY_NAMES = ("label",)
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = ("scale", "igsize", "x", "y", "u", "v", "w")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("scale", "igsize", "X", "Y", "U", "V", "W")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("scale", "igsize", "x", "y", "u", "v", "w")
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

    PREFIX = "phase"

    def __init__(self, **kwargs) -> NoReturn:
        super(Phase, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"igsize": 0., "x": 0., "y": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class PhaseL(LoopN):
    """Phases of the sample.

    Attributes
    ----------
        - label (mandatory)
        - scale, igsize, u, v, w, x, y (optional)
    """

    ITEM_CLASS = Phase
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None) -> NoReturn:
        super(PhaseL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#  loop_
#  _phase_label
#  _phase_scale
#  _phase_igsize
#   Fe3O4 0.02381() 0.2
#   Al    2.0 -0.3()
#   """

# obj = PhaseL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["Al"], end="\n\n")
