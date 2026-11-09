"""ParentSpaceGroup class (mCIF category parent_space_group)."""
from typing import NoReturn

from cryspy.B_parent_classes.cl_1_item import ItemN


class ParentSpaceGroup(ItemN):
    """Paramagnetic (parent) space group of a magnetic structure.

    The magnetic structure with a propagation vector is described by the
    operations of the group of k (``space_group_symop_magn_k``); the parent
    group is kept for information and for the magnetic domains
    (``Crystal.report_magnetic_domains``, ``Diffrn.define_magnetic_domains``).
    It is not used to calculate structure factors. Tags as in mCIF.
    """

    ATTR_MANDATORY_NAMES = ("name_hm_alt", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("name_H-M_alt", )

    ATTR_OPTIONAL_NAMES = ("it_number", "it_coordinate_system_code")
    ATTR_OPTIONAL_TYPES = (int, str)
    ATTR_OPTIONAL_CIF = ("IT_number", "IT_coordinate_system_code")

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

    PREFIX = "parent_space_group"

    def __init__(self, **kwargs) -> NoReturn:
        super(ParentSpaceGroup, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def take_space_group(self):
        """SpaceGroup object of the parent group."""
        from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
        kwargs = {"name_hm_alt": self.name_hm_alt}
        if self.is_attribute("it_coordinate_system_code"):
            kwargs["it_coordinate_system_code"] = self.it_coordinate_system_code
        space_group = SpaceGroup(**kwargs)
        space_group.form_object()
        return space_group
