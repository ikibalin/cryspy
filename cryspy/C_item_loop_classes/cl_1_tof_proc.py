import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class TOFProc(ItemN):
    """TOFProc classes.

    Attributes
    ----------
        - time_corrected, d_spacing, intensity_up_net,
          intensity_down_net, intensity_up_total, intensity_down_total, 
          intensity_bkg_calc, intensity_net, intensity_total,
          intensity_diff_total, intensity_up, intensity_up_sigma,
          intensity_down, intensity_down_sigma, intensity,
          intensity_sigma (optional)

    time is given in microseconds.
    """
    ATTR_MANDATORY_NAMES = ("time", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("time", )

    ATTR_OPTIONAL_NAMES = (
        "time_corrected", "d_spacing", "intensity_up_net",
        "intensity_down_net", "intensity_up_total", "intensity_down_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_up", "intensity_up_sigma",
        "intensity_down", "intensity_down_sigma", "intensity",
        "intensity_sigma")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float,
                           float, float, float, float, float, float, float, 
                           float, float)
    ATTR_OPTIONAL_CIF = (
        "time_corrected", "d_spacing", "intensity_up_net",
        "intensity_down_net", "intensity_up_total", "intensity_down_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_up", "intensity_up_sigma",
        "intensity_down", "intensity_down_sigma", "intensity",
        "intensity_sigma")

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

    # formats if cif format
    D_FORMATS = {
        "time_corrected": "{:.2f}", "d_spacing": "{:.5f}",
        "intensity_up_net": "{:.2f}", "intensity_down_net": "{:.2f}",
        "intensity_up_total": "{:.2f}", "intensity_down_total": "{:.2f}",
        "intensity_bkg_calc": "{:.2f}", "intensity_net": "{:.2f}",
        "intensity_total": "{:.2f}", "intensity_diff_total": "{:.2f}",
        "intensity_down": "{:.2f}", "intensity_down_sigma": "{:.2f}",
        "intensity": "{:.2f}", "intensity_sigma": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_proc"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFProc, self).__init__()

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


class TOFProcL(LoopN):
    """TOFProcL class.
    """
    ITEM_CLASS = TOFProc
    ATTR_INDEX = "time"
    def __init__(self, loop_name = None) -> NoReturn:
        super(TOFProcL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
#   loop_
#   _tof_proc_time
#   _tof_proc_time_corrected
#   _tof_proc_d_spacing
#   _tof_proc_intensity_up_net
#   _tof_proc_intensity_down_net
#   _tof_proc_intensity_up_total
#   _tof_proc_intensity_down_total
#   _tof_proc_intensity_bkg_calc
#   _tof_proc_intensity_up
#   _tof_proc_intensity_up_sigma
#   _tof_proc_intensity_down
#   _tof_proc_intensity_down_sigma
#   4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000"""

# obj = TOFProcL.from_cif(s_cont)
# print(obj, end="\n\n")
