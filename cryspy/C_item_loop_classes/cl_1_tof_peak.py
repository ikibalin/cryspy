import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFPeak(ItemN):
    """Integrated intensities for time-of-flight experiments.

    This section contains peak information extracted from the
    measured or, if present, the processed diffractogram. Each
    peak in this table will have a unique label (see _pd_peak_id).
    The reflections and phases associated with each peak will be
    specified in other sections (see the _pd_refln_ and
    _pd_phase_ sections).

    Note that peak positions are customarily determined from the
    processed diffractogram and thus corrections for position
    and intensity will have been previously applied.

    Attributes:
        - index_h, index_k, index_l (mandatory)
        - index_mult, time, intensity_plus, intensity_minus, width_time
    """
    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = ("index_mult", "time", "intensity_plus",
                           "intensity_minus", "width_time")
    ATTR_OPTIONAL_TYPES = (int, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("index_mult", "time", "intensity_plus",
                         "intensity_minus", "width_time")

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
    D_FORMATS = {"time": "{:.2f}", "intensity_plus": "{:.2f}",
                 "intensity_minus": "{:.2f}", "width_time": "{:.5f}"}

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

    PREFIX = "tof_peak"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFPeak, self).__init__()

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


class TOFPeakL(LoopN):
    """Integrated intensities for time-of-flight experiments.

    This section contains peak information extracted from the
    measured or, if present, the processed diffractogram. Each
    peak in this table will have a unique label (see _pd_peak_id).
    The reflections and phases associated with each peak will be
    specified in other sections (see the _tof_refln_ and
    _phase_ sections).

    Note that peak positions are customarily determined from the
    processed diffractogram and thus corrections for position
    and intensity will have been previously applied.

    Attributes:
        - index_h, index_k, index_l (mandatory)
        - index_mult, time, intensity_plus, intensity_minus, width_time

    """
    ITEM_CLASS = TOFPeak
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFPeakL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


# s_cont = """
#   loop_
#   _tof_peak_index_h
#   _tof_peak_index_k
#   _tof_peak_index_l
#   _tof_peak_index_mult
#   _tof_peak_time
#   _tof_peak_intensity_plus
#   _tof_peak_intensity_minus
#   _tof_peak_width_time
#   2  2  0  4 17.2 100.0  90.0  2.3
# """

# obj = TOFPeakL.from_cif(s_cont)
# print(obj, end="\n\n")
