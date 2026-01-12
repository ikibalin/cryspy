from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_3_extinction import \
    calc_extinction_2

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class JanaTwin(ItemN):
    """
    Jana Twin class.

    Attributes
    ----------
        - mosaicity 100.0 # in minutes
        - radius    50.0 # in micrometers
        - model     "gauss" or "lorentz"



    """

    ATTR_MANDATORY_NAMES = ( "id","matrix_11", "matrix_12", "matrix_13",
                            "matrix_21", "matrix_22", "matrix_23",
                            "matrix_31", "matrix_32", "matrix_33",
                            "fraction",)
    ATTR_MANDATORY_TYPES = (int, float, float, float, float, float, float, float, float, float, float)
    ATTR_MANDATORY_CIF = ("id","matrix[1][1]", "matrix[1][2]", "matrix[1][3]",
                            "matrix[2][1]", "matrix[2][2]", "matrix[2][3]",
                            "matrix[3][1]", "matrix[3][2]", "matrix[3][3]",
                            "fraction",)

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("fraction", )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

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

    PREFIX = "jana_twin"

    def __init__(self, **kwargs) -> NoReturn:
        super(JanaTwin, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"fraction": 0., }

        # defined for ani integer and float parameters
        D_MAX = {"fraction": 1., }

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class JanaTwinL(LoopN):
    """
    Description of extinction in loop.

    """
    ITEM_CLASS = JanaTwin
    ATTR_INDEX = "id"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(JanaTwinL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        res = {}
        l_fraction = self.fraction
        l_fraction_refinement = self.fraction_refinement
        res["twin_fraction"] = numpy.array(l_fraction, dtype=float)
        res["flags_twin_fraction"] = numpy.array(l_fraction_refinement, dtype=bool)
        res["twin_matrices"] = numpy.stack([
            numpy.array(self.matrix_11, dtype=float),
            numpy.array(self.matrix_12, dtype=float),
            numpy.array(self.matrix_13, dtype=float),
            numpy.array(self.matrix_21, dtype=float),
            numpy.array(self.matrix_22, dtype=float),
            numpy.array(self.matrix_23, dtype=float),
            numpy.array(self.matrix_31, dtype=float),
            numpy.array(self.matrix_32, dtype=float),
            numpy.array(self.matrix_33, dtype=float),
            ],axis=0)
        return res
