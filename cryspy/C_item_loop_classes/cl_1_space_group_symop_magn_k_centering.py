"""SpaceGroupSymopMagnKCentering, SpaceGroupSymopMagnKCenteringL classes."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class SpaceGroupSymopMagnKCentering(ItemN):
    """Lattice centring translation used with space_group_symop_magn_k.

xyz is a translation such as "x+1/2,y+1/2,z+1/2". It is a lattice
translation, so the Fourier component of the moment of the translated atom
is multiplied by exp(-2 pi i k.t) automatically. The identity "x,y,z" may
be omitted.
    """

    ATTR_MANDATORY_NAMES = ("id", "xyz")
    ATTR_MANDATORY_TYPES = (str, str)
    ATTR_MANDATORY_CIF = ("id", "xyz")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ()
    ATTR_SIGMA = ()
    ATTR_CONSTR_FLAG = ()
    ATTR_REF_FLAG = ()
    ATTR_CONSTR_MARK = ()

    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "space_group_symop_magn_k_centering"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnKCentering, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class SpaceGroupSymopMagnKCenteringL(LoopN):
    """Lattice centring translations used with space_group_symop_magn_k."""

    ITEM_CLASS = SpaceGroupSymopMagnKCentering
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnKCenteringL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
