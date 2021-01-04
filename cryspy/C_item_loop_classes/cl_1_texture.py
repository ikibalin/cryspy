from typing import NoReturn
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

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

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
    def __init__(self, loop_name = None) -> NoReturn:
        super(TextureL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#   loop_
#  _texture_g_1 
#  _texture_g_2 
#  _texture_h_ax 
#  _texture_k_ax 
#  _texture_l_ax 
#  0.1239 0.94211 -0.66119 -0.0541 3.0613
#  0.1239(27) 0.94211(18) -0.66119(41) -0.0541(6) 3.0613(7)
#  """

# obj = TextureL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj[0], end="\n\n")
