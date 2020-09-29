from typing import NoReturn, Tuple, List
import copy
import numpy
from fractions import Fraction
from cryspy.A_functions_base.function_1_strings import \
    transform_string_to_r_b, transform_r_b_to_string
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class SpaceGroupSymopMagnCentering(ItemN):
    """Centering or anti-centering translation.

    Given in a BNS-supercell description of a magnetic space group.

    Attributes
    ----------
    - xyz (mandatory)
    - description, id (optioanl)
    """

    ATTR_MANDATORY_NAMES = ("xyz", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("xyz", )

    ATTR_OPTIONAL_NAMES = ("id", "description")
    ATTR_OPTIONAL_TYPES = (str, str)
    ATTR_OPTIONAL_CIF = ("id", "description")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("r", "b", "r_11", "r_12", "r_13", "r_21", "r_22",
                      "r_23", "r_31", "r_32", "r_33", "b_1", "b_2", "b_3",
                      "theta")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "space_group_symop_magn_centering"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroupSymopMagnCentering, self).__init__()

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

    def form_object(self) -> NoReturn:
        """Form object."""
        xyz = self.xyz
        if xyz is None:
            return False
        r, b = transform_string_to_r_b(xyz, labels=("x", "y", "z"))

        self.__dict__["r_11"] = r[0, 0]
        self.__dict__["r_12"] = r[0, 1]
        self.__dict__["r_13"] = r[0, 2]
        self.__dict__["r_21"] = r[1, 0]
        self.__dict__["r_22"] = r[1, 1]
        self.__dict__["r_23"] = r[1, 2]
        self.__dict__["r_31"] = r[2, 0]
        self.__dict__["r_32"] = r[2, 1]
        self.__dict__["r_33"] = r[2, 2]
        self.__dict__["b_1"] = b[0]
        self.__dict__["b_2"] = b[1]
        self.__dict__["b_3"] = b[2]
        self.__dict__["theta"] = int(b[3])


class SpaceGroupSymopMagnCenteringL(LoopN):
    """This loop provides a list of centering or anti-centering
     translation in a BNS-supercell description of a magnetic space
     group.

     Keeping the centering and anti-centering translations in a
     separate loop leaves only representative point operations in the
     main SPACE_GROUP_SYMOP_MAGN_OPERATION loop. The direct sum of
     the two loops produces the full set of representative operations
     of the magnetic space group.  This centering loop is optional, so
     that it is always possible to include all of the symmetry
     operations in the main loop.

     When this centering loop is employed, the representative point
     operations in the main SPACE_GROUP_SYMOP_MAGN_OPERATION loop may
     not form a closed subgroup, but instead  generate some of the
     fractional translations of the centering loop.  Despite this
     annoyance, a separate centering loop is important because
     magnetic structures tend to have a relatively large number of
     centering and anti-centering translations, which can make the
     resulting list of operators very long and unintuitive, especially
     when working  in non-standard settings.

     One could argue that anti-centering operations belong in the main
     representative- point-operation loop since they are not actually
     translations of the magnetic lattice.   In fact, a pure time
     reversal is a generator of the magnetic point group of a  type-4
     magnetic space group.  Nevertheless, this centering loop is
     defined to  include the anti-centerings due to the common
     practice of referring to a "black and white" lattice of
     centerings and anti-centerings."""

    ITEM_CLASS = SpaceGroupSymopMagnCentering
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None) -> NoReturn:
        super(SpaceGroupSymopMagnCenteringL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

# s_cont = """
# loop_
# _space_group_symop_magn_centering_id
# _space_group_symop_magn_centering_xyz
# 1   x,y,z,+1
# """

# obj = SpaceGroupSymopMagnCenteringL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["1"], end="\n\n")
