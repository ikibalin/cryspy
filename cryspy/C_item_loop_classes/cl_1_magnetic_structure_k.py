"""AtomSiteMoment, AtomSiteMomentL classes are given."""
import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.A_functions_base.matrix_operations import calc_m_v

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.A_functions_base.unit_cell import calc_m_m_norm_by_unit_cell_parameters

class MagneticStructureK(ItemN):
    """MagneticStructureK class.


    Attributes
    ----------
        - number, vector_1, vector_2, vector_3 (mandatory)
    """

    ATTR_MANDATORY_NAMES = ("number", "vector_1", "vector_2", "vector_3")
    ATTR_MANDATORY_TYPES = (int, float, float, float)
    ATTR_MANDATORY_CIF = ("number", "vector_1", "vector_2", "vector_3")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

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

    PREFIX = "magnetic_structure_k"

    def __init__(self, **kwargs) -> NoReturn:
        super(MagneticStructureK, self).__init__()

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

    def to_cif(self, separator: str = ".") -> str:
        return super(MagneticStructureK, self).to_cif(separator=separator)

class MagneticStructureKL(LoopN):
    """MagneticStructureKL class.
    """

    ITEM_CLASS = MagneticStructureK
    ATTR_INDEX = "number"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MagneticStructureKL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def to_cif(self, separator: str = ".") -> str:
        return super(MagneticStructureKL, self).to_cif(separator=separator)
    

    def get_dictionary(self):
        """Get dictionary of the class attributes."""
        res = {}
        if len(self.items) > 0:
            res["magnetic_structure_k_number"] = numpy.array(self.number, dtype=int)
            res["magnetic_structure_k_vector"] = numpy.array(
                [self.vector_1, self.vector_2, self.vector_3], dtype=float)
        return res