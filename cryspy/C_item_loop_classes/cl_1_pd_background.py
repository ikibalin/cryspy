import numpy
from typing import NoReturn


from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class PdBackground(ItemN):
    """Background point in powder 1d experiment.
    """
    ATTR_MANDATORY_NAMES = ("ttheta", "intensity")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("2theta", "intensity")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("intensity", )
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

    PREFIX = "pd_background"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdBackground, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"ttheta": 0}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class PdBackgroundL(LoopN):
    """Background points in powder 1d experiment.
    """
    ITEM_CLASS = PdBackground
    ATTR_INDEX = "ttheta"
    def __init__(self, loop_name = None) -> NoReturn:
        super(PdBackgroundL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   
    def interpolate_by_points(self, tth):
        l_ttheta = self.ttheta
        l_intensity = self.intensity
        tth_b = numpy.array(l_ttheta, dtype=float)
        int_b = numpy.array(l_intensity, dtype=float)
        if len(l_ttheta) == 0:
            int_1d = numpy.zeros(tth.size, dtype=float)
        else:
            int_1d = numpy.interp(tth, tth_b, int_b)
        return int_1d

# s_cont = """
#  loop_
#  _pd_background_ttheta
#  _pd_background_intensity
#   4.5  256.0
#   40.0  158(2)
#   80.0  65.0 
# """

# obj = PdBackgroundL.from_cif(s_cont)
# print(obj, end="\n\n")
