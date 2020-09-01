import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class PdMeas(ItemN):
    """
    PdInstrReflexAsymmetry describes asymmetry of Bragg reflections for
    1d powder diffractometer.

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

    ATTR_OPTIONAL_NAMES = ("intensity_up", "intensity_up_sigma", 
                           "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("intensity_up", "intensity_up_sigma", 
                           "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")

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

    PREFIX = "pd_meas"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdMeas, self).__init__()

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


class PdMeasL(LoopN):
    """
    This section contains the measured diffractogram and information
    about the conditions used for the measurement of the diffraction 
    data set, prior to processing and application of correction
    terms. While additional information may be added to the CIF
    as data are processed and transported between laboratories
    (possibly with the addition of a new _pd_block_id entry), the
    information in this section of the CIF will rarely be changed
    once data collection is complete.

    Where possible, measurements in this section should have no
    post-collection processing applied (normalizations, corrections,
    smoothing, zero-offset corrections etc.). Such corrected
    measurements should be recorded in the _pd_proc_ section.

    Data sets that are measured as counts, where a standard
    uncertainty can be considered equivalent to the standard
    deviation and where the standard deviation can be estimated
    as the square root of the number of counts recorded, should
    use the _pd_meas_counts_ fields. All other intensity values
    should be recorded using _pd_meas_intensity_.

    """
    ITEM_CLASS = PdMeas
    ATTR_INDEX = "ttheta"
    def __init__(self, loop_name = None) -> NoReturn:
        super(PdMeasL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   
    def is_polarized(self):
        """
        True if polaraized data are defined

        and

        False for unpolaraized data
        """
        flag = False
        if len(self.items) != 0:
            if ((self.items[0].intensity_up is not None) &
                (self.items[0].intensity_down is not None)):
                flag = True
        return flag


# s_cont = """
#  loop_
#  _pd_meas_ttheta
#  _pd_meas_intensity_up
#  _pd_meas_intensity_up_sigma
#  _pd_meas_intensity_down
#  _pd_meas_intensity_down_sigma
#   4.00   465.80000   128.97000   301.88000   129.30000
#   4.20   323.78000   118.22000   206.06000   120.00000
# """

# obj = PdMeasL.from_cif(s_cont)
# print(obj, end="\n\n")
