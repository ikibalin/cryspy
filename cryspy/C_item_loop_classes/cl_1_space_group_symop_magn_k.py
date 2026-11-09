"""SpaceGroupSymopMagnK, SpaceGroupSymopMagnKL classes are given."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class SpaceGroupSymopMagnK(ItemN):
    """Symmetry operation of the group of k with its action on the moment.

Like the FullProf pair SYMM/MSYM: operation_xyz acts on the atom position
("-x+1/2,y,-z"), moment_uvw acts on the Fourier component along a, b, c
("-u,v,-w") and phase (in units of 2 pi) multiplies it by
exp(-2 pi i phase).

The operations (with the centring translations of the
space_group_symop_magn_k_centering loop) are the only symmetry of the
crystal: they generate the positions of all atoms of atom_site, nuclear
and magnetic, so every atom with a moment is an atom of the structure.
Do not give space_group together with this loop.
    """

    ATTR_MANDATORY_NAMES = ("id", "operation_xyz", "moment_uvw")
    ATTR_MANDATORY_TYPES = (str, str, str)
    ATTR_MANDATORY_CIF = ("id", "operation_xyz", "moment_uvw")

    ATTR_OPTIONAL_NAMES = ("phase", )
    ATTR_OPTIONAL_TYPES = (float, )
    ATTR_OPTIONAL_CIF = ("phase", )

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
    D_DEFAULT = {"phase": 0.}

    PREFIX = "space_group_symop_magn_k"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnK, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class SpaceGroupSymopMagnKL(LoopN):
    """Symmetry operations of the group of k with their action on moments.

See SpaceGroupSymopMagnK.
    """

    ITEM_CLASS = SpaceGroupSymopMagnK
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnKL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
