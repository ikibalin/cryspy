import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class RefineLs(ItemN):
    """
    Data items in the REFINE_LS category record details about the
    structure-refinement parameters.

    By default the UB matrix are given in ccsl notation

    Attributes
        - goodness_of_fit_all (mandatory)
        - number_reflns (mandatory)
        - number_parameters (mandatory)
        - number_restraints
        - number_constraints
        - weighting_scheme

    Internal attributes:
        - 

    Internal protected attributes:
        - 
    """
    ATTR_MANDATORY_NAMES = ("goodness_of_fit_all","number_reflns")
    ATTR_MANDATORY_TYPES = (float, int)
    ATTR_MANDATORY_CIF = ("goodness_of_fit_all", "number_reflns")

    ATTR_OPTIONAL_NAMES = ("number_restraints", "number_constraints",
                           "weighting_scheme", "number_parameters")
    ATTR_OPTIONAL_TYPES = (int, int, str, int)
    ATTR_OPTIONAL_CIF = ("number_restraints", "number_constraints",
                         "weighting_scheme", "number_parameters")

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
    D_FORMATS = {"goodness_of_fit_all": "{:.3f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"weighting_scheme": ["sigma", "unit", "calc"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "refine_ls"

    def __init__(self, **kwargs) -> NoReturn:
        super(RefineLs, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"number_reflns": 0, "goodness_of_fit_all": 0.,
                 "number_parameters": 0, "number_restraints": 0,
                 "number_constraints": 0}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)



class RefineLsL(LoopN):
    """
    Description of AtomSite in loop.

    """
    ITEM_CLASS = RefineLs
    ATTR_INDEX = "id"
    def __init__(self, loop_name = None) -> NoReturn:
        super(RefineLsL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name


# s_cont = """
#   _refine_ls_goodness_of_fit_all        2.74
#   _refine_ls_weighting_scheme           sigma
#   _refine_ls_number_reflns           1408
#   _refine_ls_number_parameters       272
#   _refine_ls_number_restraints       0
#   _refine_ls_number_constraints      0
#   """

# obj = RefineLs.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.number_reflns, end="\n\n")
