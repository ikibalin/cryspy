import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class PdProc(ItemN):
    """
    This section contains the diffraction data set after processing
    and application of correction terms. If the data set is
    reprocessed, this section may be replaced (with the addition of
    a new _pd_block_id entry).

    Mandatory attributes:
        - ub_11, ub_12, ub_13, ub_21, ub_22, ub_23, ub_31, ub_32, ub_33

    Optional attributes:
        - occupancy
        - adp_type
        - u_iso_or_equiv
        - u_equiv_geom_mean
        - b_iso_or_equiv
        - multiplicity
        - wyckoff_symbol
        - cartn_x
        - cartn_y
        - cartn_z

    Internal attributes:
        - scat_length_neutron

    Internal protected attributes:
        - space_group_wyckoff
        - constr_number
    """
    ATTR_MANDATORY_NAMES = ("ttheta", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("2theta", )

    ATTR_OPTIONAL_NAMES = (
        "ttheta_corrected", "d_spacing", "intensity_up_net",
        "intensity_down_net", "intensity_up_total", "intensity_down_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_up", "intensity_up_sigma",
        "intensity_down", "intensity_down_sigma", "intensity",
        "intensity_sigma")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float,
                           float, float, float, float, float, float, float, 
                           float, float)
    ATTR_OPTIONAL_CIF = (
        "2theta_corrected", "d_spacing", "intensity_up_net",
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
        "ttheta_corrected": "{:.2f}", "d_spacing": "{:.5f}",
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

    PREFIX = "pd_proc"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdProc, self).__init__()

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


class PdProcL(LoopN):
    """
    This section contains the diffraction data set after processing
    and application of correction terms. If the data set is
    reprocessed, this section may be replaced (with the addition of
    a new _pd_block_id entry).

    """
    ITEM_CLASS = PdProc
    ATTR_INDEX = "ttheta"
    def __init__(self, loop_name = None) -> NoReturn:
        super(PdProcL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
#  loop_
#  _pd_proc_ttheta
#  _pd_proc_ttheta_corrected
#  _pd_proc_d_spacing
#  _pd_proc_intensity_up_net
#  _pd_proc_intensity_down_net
#  _pd_proc_intensity_up_total
#  _pd_proc_intensity_down_total
#  _pd_proc_intensity_bkg_calc
#  _pd_proc_intensity_up
#  _pd_proc_intensity_up_sigma
#  _pd_proc_intensity_down
#  _pd_proc_intensity_down_sigma
#   4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000"""

# obj = PdProcL.from_cif(s_cont)
# print(obj, end="\n\n")
