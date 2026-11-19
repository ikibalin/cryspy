"""GeomBond, GeomAngle, GeomTorsion and their loops (core CIF categories
geom_bond, geom_angle, geom_torsion): interatomic distances, bond angles
and torsion angles of a crystal.

Site symmetry codes 'n_klm' (core CIF dictionary): n is the number of the
symmetry operation (centring included) as Crystal.report_distances_and_angles lists
it in the loop _space_group_symop, k, l, m are 5 + the lattice translation
along a, b, c; '.' is the atom of atom_site itself.
"""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class GeomBond(ItemN):
    """Interatomic distance (Å) between two atoms of atom_site.

The value is written with its standard uncertainty, '1.4568(12)', as in a
publication CIF. It is calculated from the structure by
``Crystal.add_distances_and_angles`` (again after a refinement), it is not a parameter
of the refinement.
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2",)
    ATTR_MANDATORY_TYPES = (str, str,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2",)

    ATTR_OPTIONAL_NAMES = ("distance", "site_symmetry_1", "site_symmetry_2", "publ_flag",)
    ATTR_OPTIONAL_TYPES = (str, str, str, str,)
    ATTR_OPTIONAL_CIF = ("distance", "site_symmetry_1", "site_symmetry_2", "publ_flag",)

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

    D_FORMATS = {}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "geom_bond"

    def __init__(self, **kwargs) -> NoReturn:
        super(GeomBond, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class GeomBondL(LoopN):
    """Interatomic distances of the crystal."""

    ITEM_CLASS = GeomBond
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(GeomBondL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class GeomAngle(ItemN):
    """Bond angle (degrees) at the second atom (the vertex).

The value is written with its standard uncertainty, '109.04(9)', as in a
publication CIF. It is calculated from the structure by
``Crystal.add_distances_and_angles`` (again after a refinement), it is not a parameter
of the refinement.
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3",)
    ATTR_MANDATORY_TYPES = (str, str, str,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3",)

    ATTR_OPTIONAL_NAMES = ("value", "site_symmetry_1", "site_symmetry_2", "site_symmetry_3",
                           "publ_flag",)
    ATTR_OPTIONAL_TYPES = (str, str, str, str, str,)
    # the value is _geom_angle itself in core CIF (DDL1)
    ATTR_OPTIONAL_CIF = ("", "site_symmetry_1", "site_symmetry_2", "site_symmetry_3",
                         "publ_flag",)

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

    D_FORMATS = {}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "geom_angle"

    def __init__(self, **kwargs) -> NoReturn:
        super(GeomAngle, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class GeomAngleL(LoopN):
    """Bond angles of the crystal."""

    ITEM_CLASS = GeomAngle
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(GeomAngleL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class GeomTorsion(ItemN):
    """Torsion angle (degrees) of four atoms, IUPAC sign.

The value is written with its standard uncertainty, '-61.2(3)', as in a
publication CIF. It is calculated from the structure by
``Crystal.add_distances_and_angles`` (again after a refinement), it is not a parameter
of the refinement.
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "atom_site_label_4",)
    ATTR_MANDATORY_TYPES = (str, str, str, str,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "atom_site_label_4",)

    ATTR_OPTIONAL_NAMES = ("value", "site_symmetry_1", "site_symmetry_2", "site_symmetry_3",
                           "site_symmetry_4", "publ_flag",)
    ATTR_OPTIONAL_TYPES = (str, str, str, str, str, str,)
    # the value is _geom_torsion itself in core CIF (DDL1)
    ATTR_OPTIONAL_CIF = ("", "site_symmetry_1", "site_symmetry_2", "site_symmetry_3",
                         "site_symmetry_4", "publ_flag",)

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

    D_FORMATS = {}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "geom_torsion"

    def __init__(self, **kwargs) -> NoReturn:
        super(GeomTorsion, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class GeomTorsionL(LoopN):
    """Torsion angles of the crystal."""

    ITEM_CLASS = GeomTorsion
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(GeomTorsionL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class GeomHbond(ItemN):
    """Hydrogen bond D–H···A: donor D, hydrogen (or deuterium) H, acceptor A.

Distances D–H, H···A, D···A in Å and the angle D–H···A in degrees, written
with their standard uncertainties. Calculated from the structure by
``Crystal.add_distances_and_angles``, not parameters of the refinement.
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_d", "atom_site_label_h", "atom_site_label_a",)
    ATTR_MANDATORY_TYPES = (str, str, str,)
    ATTR_MANDATORY_CIF = ("atom_site_label_d", "atom_site_label_h", "atom_site_label_a",)

    ATTR_OPTIONAL_NAMES = ("distance_dh", "distance_ha", "distance_da", "angle_dha",
                           "site_symmetry_d", "site_symmetry_h", "site_symmetry_a",
                           "publ_flag",)
    ATTR_OPTIONAL_TYPES = (str, str, str, str, str, str, str, str,)
    # CIF is case-insensitive: _geom_hbond_atom_site_label_D is read too
    ATTR_OPTIONAL_CIF = ("distance_dh", "distance_ha", "distance_da", "angle_dha",
                         "site_symmetry_d", "site_symmetry_h", "site_symmetry_a",
                         "publ_flag",)

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

    D_FORMATS = {}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "geom_hbond"

    def __init__(self, **kwargs) -> NoReturn:
        super(GeomHbond, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class GeomHbondL(LoopN):
    """Hydrogen bonds of the crystal."""

    ITEM_CLASS = GeomHbond
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(GeomHbondL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
