from typing import NoReturn

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


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
        "chi_32_calc", "chi_33_calc", "moment_11_calc", "moment_12_calc",
        "moment_13_calc", "moment_21_calc", "moment_22_calc", "moment_23_calc", 
        "moment_31_calc", "moment_32_calc", "moment_33_calc")
    ATTR_OPTIONAL_TYPES = (float, float, complex, complex, complex,
                           complex, complex, complex, complex, complex, 
                           complex, complex, complex, complex, complex,
                           complex, complex, complex, complex, complex)
    ATTR_OPTIONAL_CIF = (
        "sintlambda", "d_spacing", "chi_11_calc", "chi_12_calc", "chi_13_calc", 
        "chi_21_calc", "chi_22_calc", "chi_23_calc", "chi_31_calc",
        "chi_32_calc", "chi_33_calc", "moment_11_calc", "moment_12_calc",
        "moment_13_calc", "moment_21_calc", "moment_22_calc", "moment_23_calc", 
        "moment_31_calc", "moment_32_calc", "moment_33_calc")

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
    D_FORMATS = {"sintlambda": "{:.5f}", "d_spacing": "{:.5f}",
                 "chi_11_calc": "{:.2f}", "chi_12_calc": "{:.2f}",
                 "chi_13_calc": "{:.2f}", "chi_21_calc": "{:.2f}",
                 "chi_22_calc": "{:.2f}", "chi_23_calc": "{:.2f}",
                 "chi_31_calc": "{:.2f}", "chi_32_calc": "{:.2f}",
                 "chi_33_calc": "{:.2f}", "moment_11_calc": "{:.2f}",
                 "moment_12_calc": "{:.2f}", "moment_13_calc": "{:.2f}",
                 "moment_21_calc": "{:.2f}", "moment_22_calc": "{:.2f}",
                 "moment_23_calc": "{:.2f}", "moment_31_calc": "{:.2f}",
                 "moment_32_calc": "{:.2f}", "moment_33_calc": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

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
    def __init__(self, loop_name = None) -> NoReturn:
        super(ReflnSusceptibilityL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name


# s_cont = """
#  loop_
#  _refln_susceptibility.index_h    
#  _refln_susceptibility.index_k    
#  _refln_susceptibility.index_l    
#  _refln_susceptibility.d_spacing       
#  _refln_susceptibility.sintlambda
#  _refln_susceptibility.chi_11_calc     
#  _refln_susceptibility.chi_12_calc     
#  _refln_susceptibility.chi_13_calc     
#  _refln_susceptibility.chi_21_calc     
#  _refln_susceptibility.chi_22_calc     
#  _refln_susceptibility.chi_23_calc     
#  _refln_susceptibility.chi_31_calc     
#  _refln_susceptibility.chi_32_calc     
#  _refln_susceptibility.chi_33_calc     
#  _refln_susceptibility.moment_11_calc  
#  _refln_susceptibility.moment_12_calc  
#  _refln_susceptibility.moment_13_calc  
#  _refln_susceptibility.moment_21_calc  
#  _refln_susceptibility.moment_22_calc  
#  _refln_susceptibility.moment_23_calc  
#  _refln_susceptibility.moment_31_calc  
#  _refln_susceptibility.moment_32_calc  
#  _refln_susceptibility.moment_33_calc  
#  2 0 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j
#  0 2 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j
#   """

# obj = ReflnSusceptibilityL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.sintlambda, end="\n\n")
# print(obj.numpy_chi_33_calc, end="\n\n")
