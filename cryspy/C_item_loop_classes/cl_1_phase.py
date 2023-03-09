"""Description of Phase and PhaseL classes."""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Phase(ItemN):
    """Phase class.

    Attributes
    ----------
        - label (mandatory)
        - scale, igsize, u, v, w, x, y (optional)
    """

    ATTR_MANDATORY_NAMES = ("label",)
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = ("scale", "igsize", "x", "y", "u", "v", "w")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("scale", "igsize", "X", "Y", "U", "V", "W")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("scale", "igsize", "x", "y", "u", "v", "w")
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

    PREFIX = "phase"

    def __init__(self, **kwargs) -> NoReturn:
        super(Phase, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"igsize": 0., "x": 0., "y": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        res = {}
        res["phase_name"] = numpy.array([self.label,], dtype=str)
        n_items = 1
        if self.is_attribute("u"):
            p_u = numpy.array([self.u,], dtype=float)
            r_u = numpy.array([self.u_refinement,], dtype=bool)
        else:
            p_u = numpy.zeros((n_items,), dtype=float)
            r_u = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("v"):
            p_v = numpy.array([self.v,], dtype=float)
            r_v = numpy.array([self.v_refinement,], dtype=bool)
        else:
            p_v = numpy.zeros((n_items,), dtype=float)
            r_v = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("w"):
            p_w = numpy.array([self.w,], dtype=float)
            r_w = numpy.array([self.w_refinement,], dtype=bool)
        else:
            p_w = numpy.zeros((n_items,), dtype=float)
            r_w = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("x"):
            p_x = numpy.array([self.x,], dtype=float)
            r_x = numpy.array([self.x_refinement,], dtype=bool)
        else:
            p_x = numpy.zeros((n_items,), dtype=float)
            r_x = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("y"):
            p_y = numpy.array([self.y,], dtype=float)
            r_y = numpy.array([self.y_refinement,], dtype=bool)
        else:
            p_y = numpy.zeros((n_items,), dtype=float)
            r_y = numpy.zeros((n_items,), dtype=bool)

        res["phase_resolution_parameters"] = numpy.stack([p_u, p_v, p_w, p_x, p_y], axis=0)
        res["flags_phase_resolution_parameters"] = numpy.stack([r_u, r_v, r_w, r_x, r_y], axis=0)

        if self.is_attribute("igsize"):
            res["phase_ig"] = numpy.array([self.igsize,], dtype=float)
            res["flags_phase_ig"] = numpy.array([self.igsize_refinement,], dtype=bool)
        else:
            res["phase_ig"] = numpy.zeros((n_items,), dtype=float)
            res["flags_phase_ig"] = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("scale"):
            res["phase_scale"] = numpy.array([self.scale,], dtype=float)
            res["flags_phase_scale"] = numpy.array([self.scale_refinement,], dtype=bool)
        else:
            res["phase_scale"] = numpy.zeros((n_items,), dtype=float)
            res["flags_phase_scale"] = numpy.zeros((n_items,), dtype=bool)
        return res


class PhaseL(LoopN):
    """Phases of the sample.

    Attributes
    ----------
        - label (mandatory)
        - scale, igsize, u, v, w, x, y (optional)
    """

    ITEM_CLASS = Phase
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(PhaseL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        res = {}
        res["phase_name"] = numpy.array(self.label, dtype=str)
        n_items = len(self.items)
        if self.is_attribute("u"):
            p_u = numpy.array(self.u, dtype=float)
            r_u = numpy.array(self.u_refinement, dtype=bool)
        else:
            p_u = numpy.zeros((n_items,), dtype=float)
            r_u = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("v"):
            p_v = numpy.array(self.v, dtype=float)
            r_v = numpy.array(self.v_refinement, dtype=bool)
        else:
            p_v = numpy.zeros((n_items,), dtype=float)
            r_v = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("w"):
            p_w = numpy.array(self.w, dtype=float)
            r_w = numpy.array(self.w_refinement, dtype=bool)
        else:
            p_w = numpy.zeros((n_items,), dtype=float)
            r_w = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("x"):
            p_x = numpy.array(self.x, dtype=float)
            r_x = numpy.array(self.x_refinement, dtype=bool)
        else:
            p_x = numpy.zeros((n_items,), dtype=float)
            r_x = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("y"):
            p_y = numpy.array(self.y, dtype=float)
            r_y = numpy.array(self.y_refinement, dtype=bool)
        else:
            p_y = numpy.zeros((n_items,), dtype=float)
            r_y = numpy.zeros((n_items,), dtype=bool)

        res["phase_resolution_parameters"] = numpy.stack([p_u, p_v, p_w, p_x, p_y], axis=0)
        res["flags_phase_resolution_parameters"] = numpy.stack([r_u, r_v, r_w, r_x, r_y], axis=0)

        if self.is_attribute("igsize"):
            res["phase_ig"] = numpy.array(self.igsize, dtype=float)
            res["flags_phase_ig"] = numpy.array(self.igsize_refinement, dtype=bool)
        else:
            res["phase_ig"] = numpy.zeros((n_items,), dtype=float)
            res["flags_phase_ig"] = numpy.zeros((n_items,), dtype=bool)

        if self.is_attribute("scale"):
            res["phase_scale"] = numpy.array(self.scale, dtype=float)
            res["flags_phase_scale"] = numpy.array(self.scale_refinement, dtype=bool)
        else:
            res["phase_scale"] = numpy.zeros((n_items,), dtype=float)
            res["flags_phase_scale"] = numpy.zeros((n_items,), dtype=bool)
        return res
