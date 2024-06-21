import numpy
from typing import NoReturn
import scipy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeasL

na = numpy.newaxis


class TOFBackgroundPoint(ItemN):
    """Background point in powder 1d experiment.
    """
    ATTR_MANDATORY_NAMES = ("time", "intensity")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("time", "intensity")

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

    PREFIX = "tof_backgroundpoint"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFBackgroundPoint, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"time": 0}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class TOFBackgroundPointL(LoopN):
    """Background points in powder 1d experiment.
    """
    ITEM_CLASS = TOFBackgroundPoint
    ATTR_INDEX = "time"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFBackgroundPointL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(
            self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def interpolate_by_points(self, np_time):
        l_time = self.time
        l_intensity = self.intensity
        time_b = numpy.array(l_time, dtype=float)
        int_b = numpy.array(l_intensity, dtype=float)
        if len(l_time) == 0:
            int_1d = numpy.zeros(np_time.size, dtype=float)
        else:
            int_1d = numpy.interp(np_time, time_b, int_b)
        return int_1d

    def get_dictionary(self):
        res = {}
        res["background_time"] = numpy.array(self.time, dtype=float)
        res["background_intensity"] = numpy.array(self.intensity, dtype=float)
        res["flags_background_intensity"] = numpy.array(
            self.intensity_refinement, dtype=bool)
        return res
