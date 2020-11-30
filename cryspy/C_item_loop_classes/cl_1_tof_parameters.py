import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class TOFParameters(ItemN):
    """
    Reflexion positions parameters in T.O.F patterns

    Attributes
    ----------
        - zero, dtt1, tsinth (mandatory)
        - dtt2 (optional)

    TOF = zero + dtt1 d + dtt2 d2

    Value of for the detector bank. Used for obtaining the wavelengths and
    for Lorentz factor correction.
    """
    ATTR_MANDATORY_NAMES = ("zero", "dtt1", "ttheta_bank")
    ATTR_MANDATORY_TYPES = (float, float, float)
    ATTR_MANDATORY_CIF = ("Zero", "Dtt1", "2theta_bank")

    ATTR_OPTIONAL_NAMES = ("dtt2", )
    ATTR_OPTIONAL_TYPES = (float, )
    ATTR_OPTIONAL_CIF = ("Dtt2", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("zero", "dtt1", "dtt2")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"zero": "{:.5f}", "dtt1": "{:.5f}", "dtt2": "{:.5f}",
                 "ttheta_bank": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"zero": 0., "dtt2": 0.}

    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_parameters"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFParameters, self).__init__()

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



class TOFParametersL(LoopN):
    """
    Description of AtomSite in loop.

    """
    ITEM_CLASS = TOFParameters
    ATTR_INDEX = None
    def __init__(self, loop_name = None) -> NoReturn:
        super(TOFParametersL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
# _tof_parameters_zero 2.921
# _tof_parameters_dtt1 6167.247
# _tof_parameters_dtt1 -2.280
# _tof_parameters_2theta_bank 145.
# """

# obj = TOFParameters.from_cif(s_cont)
# print(obj, end="\n\n")
