"""Description of classes Refln, ReflnL."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Refln(ItemN):
    """Nuclear structure factor of reflection.

    Data items in the REFLN category record details about the
    reflections used to determine the ATOM_SITE data items.
    The REFLN data items refer to individual reflections and must
    be included in looped lists.

    The REFLNS data items specify the parameters that apply to all
    reflections. The REFLNS data items are not looped.

    Mandatory attributes:
        - index_h, index_k, index_l

    Optional attributes:
        - a_calc, a_meas, b_calc, b_meas, class_code, crystal_id,
        - d_spacing, f_calc, f_meas, f_sigma, f_squared_calc,
        - f_squared_meas, f_squared_sigma, include_status,
        - intensity_calc, intensity_meas, intensity_sigma,
        - mean_path_length_tbar, phase_calc, phase_meas,
        - refinement_status, scale_group_code, sintlambda,
        - symmetry_epsilon, symmetry_multiplicity, wavelength,
        - wavelength_id

    """

    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = (
        "a_calc", "a_meas", "b_calc", "b_meas", "class_code", "crystal_id",
        "d_spacing", "f_calc", "f_meas", "f_sigma", "f_squared_calc", "f_polarized_squared_calc",
        "f_squared_meas", "f_squared_sigma", "include_status",
        "intensity_calc", "intensity_meas", "intensity_sigma",
        "mean_path_length_tbar", "phase_calc", "phase_meas",
        "refinement_status", "scale_group_code", "sintlambda",
        "symmetry_epsilon", "symmetry_multiplicity", "wavelength",
        "wavelength_id")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, str, str,
                           float, complex, complex, complex, float, float,
                           float, float, str, float, float, float,
                           float, str, str, str, str, float,
                           str, int, float, str)
    ATTR_OPTIONAL_CIF = (
        "a_calc", "a_meas", "b_calc", "b_meas", "class_code", "crystal_id",
        "d_spacing", "f_calc", "f_meas", "f_sigma", "f_squared_calc", "f_polarized_squared_calc",
        "f_squared_meas", "f_squared_sigma", "include_status",
        "intensity_calc", "intensity_meas", "intensity_sigma",
        "mean_path_length_tbar", "phase_calc", "phase_meas",
        "refinement_status", "scale_group_code", "sint/lambda",
        "symmetry_epsilon", "symmetry_multiplicity", "wavelength",
        "wavelength_id")

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
    D_FORMATS = {"f_calc": "{:.2f}", "a_calc": "{:.4f}", "b_calc": "{:.4f}",
                 "f_squared_calc": "{:.2f}", "f_polarized_squared_calc": "{:.2f}",
                 "sintlambda": "{:.5f}", "d_spacing": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"include_status": ["o", "<", "-", "x", "h", "r"],
                     "refinement_status": ["incl", "excl", "extn"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "refln"

    def __init__(self, **kwargs) -> NoReturn:
        super(Refln, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"sintlambda": 0., "d_spacing": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class ReflnL(LoopN):
    """Nuclear structure factors for reflections."""

    ITEM_CLASS = Refln
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ReflnL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


# s_cont = """
#   loop_
#   _refln_index_h
#   _refln_index_k
#   _refln_index_l
#   _refln_d_spacing
#   _refln_A_calc
#   _refln_B_calc
#   _refln_f_calc
#   0 0 2 2.315 3.25  1.232 3.25+1.232j
#   2 2 0 4.213 5.00 -4.05 4.5
#   """

# obj = ReflnL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.numpy_f_calc, end="\n\n")

# import numpy

# index_h = numpy.array([1,2,1,2], dtype=int)
# index_k = numpy.array([2,0,2,0], dtype=int)
# index_l = numpy.array([1,1,1,0], dtype=int)

# obj_2 = ReflnL()
# obj_2.numpy_index_h = index_h
# obj_2.numpy_index_k = index_k
# obj_2.numpy_index_l = index_l
# obj_2.numpy_to_items()
# print(obj_2)
