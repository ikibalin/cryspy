from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.A_functions_base.matrix_operations import calc_m1_m2, calc_vt_m_v

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import DiffrnOrientMatrix

class ReflnSusceptibility(ItemN):
    """Susceptibility structure factor tensor.
    
    Data items in the REFLN_SUSCEPTIBILITY category record details about the
    susceptibility structure factor tensor and magnetization structure factor 
    tensor calculated for different reflections.

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
    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = (
        "sintlambda", "d_spacing", "chi_11_calc", "chi_12_calc", "chi_13_calc", 
        "chi_21_calc", "chi_22_calc", "chi_23_calc", "chi_31_calc",
        "chi_32_calc", "chi_33_calc", )
    ATTR_OPTIONAL_TYPES = (float, float, complex, complex, complex,
                           complex, complex, complex, complex, complex, 
                           complex, complex, complex, complex, )
    ATTR_OPTIONAL_CIF = (
        "sint/lambda", "d_spacing", "chi_11_calc", "chi_12_calc", "chi_13_calc", 
        "chi_21_calc", "chi_22_calc", "chi_23_calc", "chi_31_calc",
        "chi_32_calc", "chi_33_calc", )

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
    D_FORMATS = {"sintlambda": "{:.5f}", "d_spacing": "{:.5f}",
                 "chi_11_calc": "{:.2f}", "chi_12_calc": "{:.2f}",
                 "chi_13_calc": "{:.2f}", "chi_21_calc": "{:.2f}",
                 "chi_22_calc": "{:.2f}", "chi_23_calc": "{:.2f}",
                 "chi_31_calc": "{:.2f}", "chi_32_calc": "{:.2f}",
                 "chi_33_calc": "{:.2f}", }

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

    PREFIX = "refln_susceptibility"

    def __init__(self, **kwargs) -> NoReturn:
        super(ReflnSusceptibility, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"sintlambda": 0, "d_spacing": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class ReflnSusceptibilityL(LoopN):
    """Susceptibility structure factor tensors.
    """
    ITEM_CLASS = ReflnSusceptibility
    ATTR_INDEX = "id"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ReflnSusceptibilityL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def calc_f_mag(self, diffrn_orient_matrix: DiffrnOrientMatrix):
        e_up = diffrn_orient_matrix.calc_e_up(phi=0, omega=0, chi=0)
        m1_ij = numpy.array([
            self.chi_11_calc, self.chi_12_calc, self.chi_13_calc,
            self.chi_21_calc, self.chi_22_calc, self.chi_23_calc,
            self.chi_31_calc, self.chi_32_calc, self.chi_33_calc], dtype=complex)
        
        m2_ij =numpy.conjugate(numpy.stack([
            m1_ij[0], m1_ij[3], m1_ij[6], m1_ij[1], m1_ij[4], m1_ij[7],
            m1_ij[2], m1_ij[5], m1_ij[8]], axis=0))
        
        m_3_ij = (calc_m1_m2(m2_ij, m1_ij))
        f_mag_sq = calc_vt_m_v(m_3_ij, e_up).real
        f_mag = numpy.sqrt(f_mag_sq)
        return f_mag

