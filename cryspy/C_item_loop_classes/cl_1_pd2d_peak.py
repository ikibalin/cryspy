"""Classes Pd2dPeak, Pd2dPeakL."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Pd2dPeak(ItemN):
    """
    Pd2dPeak class.

    This section contains peak information extracted from the
    measured or, if present, the processed diffractogram.

    Attributes
    ----------
        - index_h, index_k, index_l (mandatory)
        - index_mult, ttheta, f_nucl_sq, f_m_p_sin_sq
        - f_m_p_cos_sq, cross_sin, width_ttheta
    """

    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = ("index_mult", "ttheta", "f_nucl_sq", "f_m_p_sin_sq",
                           "f_m_p_cos_sq", "cross_sin", "width_ttheta",
                           "sintlambda", "d_spacing")
    ATTR_OPTIONAL_TYPES = (int, float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("index_mult", "2theta", "f_nucl_sq", "f_m_p_sin_sq",
                         "f_m_p_cos_sq", "cross_sin", "width_2theta",
                         "sint/lambda", "d_spacing")

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
    D_FORMATS = {"ttheta": "{:.2f}", "f_nucl_sq": "{:.2f}",
                 "f_m_p_sin_sq": "{:.2f}", "f_m_p_cos_sq": "{:.2f}",
                 "cross_sin": "{:.2f}", "width_ttheta": "{:.5f}",
                 "sintlambda": "{:.5f}", "d_spacing": "{:.5f}"}

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

    PREFIX = "pd2d_peak"

    def __init__(self, **kwargs) -> NoReturn:
        super(Pd2dPeak, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"ttheta": 0., "width_ttheta": 0.,
                 "sintlambda": 0., "d_spacing": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class Pd2dPeakL(LoopN):
    """Pd2dPeakL class."""

    ITEM_CLASS = Pd2dPeak
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(Pd2dPeakL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


# s_cont = """
#  loop_
#  _pd2d_peak_index_h
#  _pd2d_peak_index_k
#  _pd2d_peak_index_l
#  _pd2d_peak_index_mult
#  _pd2d_peak_ttheta
#  _pd2d_peak_f_nucl_sq
#  _pd2d_peak_f_m_p_sin_sq
#  _pd2d_peak_f_m_p_cos_sq
#  _pd2d_peak_cross_sin
#  _pd2d_peak_width_ttheta
#   2  2  0  4 17.2 100.0 101.2   90.0  87.4 2.3
#   2  1  3  2 19.1  78.6 101.0   92.0  82.7 5.7
# """

# obj = Pd2dPeakL.from_cif(s_cont)
# print(obj, end="\n\n")
