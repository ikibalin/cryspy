"""SpaceGroupSymopMagnSsgOperation, SpaceGroupSymopMagnSsgOperationL classes are given (mCIF category space_group_symop_magn_ssg_operation)."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class SpaceGroupSymopMagnSsgOperation(ItemN):
    """Operation of the magnetic superspace group (mCIF).

``algebraic`` is written as ``x1,x2,x3,x4,θ``, e.g. ``-x1+1/2,x2,-x3,-x4+1/2,-1``:
the action on the coordinates x1, x2, x3, on the internal coordinate x4
(``x4 = k.r`` for the moment of the atom at r) and the time reversal θ = ±1.
These operations are the only symmetry of the crystal: they generate the
nuclear and the magnetic structure together, so ``space_group`` is not
given (the paramagnetic group can be kept as ``parent_space_group``).
    """

    ATTR_MANDATORY_NAMES = ("algebraic", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("algebraic", )

    ATTR_OPTIONAL_NAMES = ("id", )
    ATTR_OPTIONAL_TYPES = (str, )
    ATTR_OPTIONAL_CIF = ("id", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {}

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "space_group_symop_magn_ssg_operation"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnSsgOperation, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class SpaceGroupSymopMagnSsgOperationL(LoopN):
    """Operations of the magnetic superspace group (mCIF).
    """

    ITEM_CLASS = SpaceGroupSymopMagnSsgOperation
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnSsgOperationL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
