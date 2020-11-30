import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class TOFBackground(ItemN):
    """TOFBackground class.

    Attributes
    ----------
        - coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8,
          coeff9, coeff10, coeff11, coeff12, coeff13, coeff14, coeff15,
          coeff16, coeff17, coeff18, id (optional)
    """
    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = (
        "coeff1", "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7",
        "coeff8", "coeff9", "coeff10", "coeff11", "coeff12", "coeff13",
        "coeff14", "coeff15", "coeff16", "coeff17", "coeff18", "id")

    ATTR_OPTIONAL_TYPES = (
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, str)
    ATTR_OPTIONAL_CIF = (
        "coeff1", "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7",
        "coeff8", "coeff9", "coeff10", "coeff11", "coeff12", "coeff13",
        "coeff14", "coeff15", "coeff16", "coeff17", "coeff18", "id")


    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = (
        "coeff1", "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7",
        "coeff8", "coeff9", "coeff10", "coeff11", "coeff12", "coeff13",
        "coeff14", "coeff15", "coeff16", "coeff17", "coeff18")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'coeff1': "{:.5f}", 'coeff2': "{:.5f}", 'coeff3': "{:.5f}",
                 'coeff4': "{:.5f}", 'coeff5': "{:.5f}", 'coeff6': "{:.5f}",
                 'coeff7': "{:.5f}", 'coeff8': "{:.5f}", 'coeff9': "{:.5f}",
                 'coeff10': "{:.5f}", 'coeff11': "{:.5f}", 'coeff12': "{:.5f}",
                 'coeff13': "{:.5f}", 'coeff14': "{:.5f}", 'coeff15': "{:.5f}",
                 'coeff16': "{:.5f}", 'coeff17': "{:.5f}", 'coeff18': "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_background"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFBackground, self).__init__()

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


class TOFBackgroundL(LoopN):
    """
    Description of TOFBackgroundL in loop.

    """
    ITEM_CLASS = TOFBackground
    ATTR_INDEX = "id"
    def __init__(self, loop_name = None) -> NoReturn:
        super(TOFBackgroundL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
# _tof_background_coeff1 24832.850
# _tof_background_coeff2 6139.244
# _tof_background_coeff3 8063.472
# _tof_background_coeff4 3125.050
# _tof_background_coeff5 2566.956
# _tof_background_coeff6 311.077
# _tof_background_coeff7 837.348
# _tof_background_coeff8 -103.742
# _tof_background_coeff9 -11.806
# """

# obj = TOFBackground.from_cif(s_cont)
# print(obj, end="\n\n")
