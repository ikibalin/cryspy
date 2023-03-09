from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Texture(ItemN):
    """
    Texture description by Modified March’s function.

    Mandatory Attributes
    --------------------
        - g_1, g_2, h_ax, k_ax, l_ax

    """
    ATTR_MANDATORY_NAMES = ("g_1", "g_2", "h_ax", "k_ax", "l_ax")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float)
    ATTR_MANDATORY_CIF = ("g_1", "g_2", "h_ax", "k_ax", "l_ax")

    ATTR_OPTIONAL_NAMES = ("label", )
    ATTR_OPTIONAL_TYPES = (str, )
    ATTR_OPTIONAL_CIF = ("label", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("g_1", "g_2", "h_ax", "k_ax", "l_ax")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'g_1': "{:.5f}", 'g_2': "{:.5f}", 'h_ax': "{:.3f}",
                 'k_ax': "{:.3f}", 'l_ax': "{:.3f}"}

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

    PREFIX = "texture"

    def __init__(self, **kwargs) -> NoReturn:
        super(Texture, self).__init__()

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


class TextureL(LoopN):
    """
    Description of Texture in loop.

    """
    ITEM_CLASS = Texture
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TextureL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        res = {}
        res["texture_name"] = numpy.array(self.label, dtype=str)
        res["texture_g1"] = numpy.array(self.g_1, dtype=float)
        res["texture_g2"] = numpy.array(self.g_2, dtype=float)
        res["texture_axis"] = numpy.array(
            [self.h_ax, self.k_ax, self.l_ax], dtype=float)
        res["flags_texture_g1"] = numpy.array(self.g_1_refinement, dtype=bool)
        res["flags_texture_g2"] = numpy.array(self.g_2_refinement, dtype=bool)
        res["flags_texture_axis"] = numpy.array(
            [self.h_ax_refinement,
                self.k_ax_refinement,
                self.l_ax_refinement], dtype=bool)
        return res