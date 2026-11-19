"""RestrDistance, RestrAngle, RestrTorsion and their loops (core CIF
categories restr_distance, restr_angle, restr_torsion): restraints of the
geometry (soft constraints, SHELX DFIX/DANG style) for the refinement.
"""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class RestrDistance(ItemN):
    """Restraint of the distance (Å) between two atoms.

``target`` and ``target_weight_param`` (its standard uncertainty sigma)
are given by the user and never recalculated; the refinement adds
((value - target) / sigma)^2 to chi_sq. ``diff`` = value - target is
calculated (Crystal.update_restraints, after every refinement). The atoms
are equivalent positions of atom_site atoms: the symmetry operation n of
the site symmetry code 'n_klm' (as in geom_bond) is kept, the lattice
translation is the one that puts the atom closest to the previous atom of
the chain (to the middle atom for an angle).
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2", "target", "target_weight_param",)
    ATTR_MANDATORY_TYPES = (str, str,  float, float,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2", "target", "target_weight_param",)

    ATTR_OPTIONAL_NAMES = ("site_symmetry_1", "site_symmetry_2", "diff", "details",)
    ATTR_OPTIONAL_TYPES = (str, str,  float, str,)
    ATTR_OPTIONAL_CIF = ("site_symmetry_1", "site_symmetry_2", "diff", "details",)

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

    D_FORMATS = {"target": "{:.4f}", "target_weight_param": "{:.4f}", "diff": "{:.4f}"}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "restr_distance"

    def __init__(self, **kwargs) -> NoReturn:
        super(RestrDistance, self).__init__()
        self.__dict__["D_MIN"] = {"target_weight_param": 0.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_labels(self) -> list:
        return [getattr(self, f"atom_site_label_{i:}") for i in range(1, 2 + 1)]

    def get_site_symmetry(self) -> list:
        # '.' (the atom itself) is read from CIF as None
        return [str(getattr(self, f"site_symmetry_{i:}", None) or ".")
                for i in range(1, 2 + 1)]


class RestrDistanceL(LoopN):
    """Distance restraints of the crystal."""

    ITEM_CLASS = RestrDistance
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(RestrDistanceL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class RestrAngle(ItemN):
    """Restraint of the bond angle (degrees) at the second atom.

``target`` and ``target_weight_param`` (its standard uncertainty sigma)
are given by the user and never recalculated; the refinement adds
((value - target) / sigma)^2 to chi_sq. ``diff`` = value - target is
calculated (Crystal.update_restraints, after every refinement). The atoms
are equivalent positions of atom_site atoms: the symmetry operation n of
the site symmetry code 'n_klm' (as in geom_bond) is kept, the lattice
translation is the one that puts the atom closest to the previous atom of
the chain (to the middle atom for an angle).
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "target", "target_weight_param",)
    ATTR_MANDATORY_TYPES = (str, str, str,  float, float,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "target", "target_weight_param",)

    ATTR_OPTIONAL_NAMES = ("site_symmetry_1", "site_symmetry_2", "site_symmetry_3", "diff", "details",)
    ATTR_OPTIONAL_TYPES = (str, str, str,  float, str,)
    ATTR_OPTIONAL_CIF = ("site_symmetry_1", "site_symmetry_2", "site_symmetry_3", "diff", "details",)

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

    D_FORMATS = {"target": "{:.2f}", "target_weight_param": "{:.2f}", "diff": "{:.2f}"}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "restr_angle"

    def __init__(self, **kwargs) -> NoReturn:
        super(RestrAngle, self).__init__()
        self.__dict__["D_MIN"] = {"target_weight_param": 0.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_labels(self) -> list:
        return [getattr(self, f"atom_site_label_{i:}") for i in range(1, 3 + 1)]

    def get_site_symmetry(self) -> list:
        # '.' (the atom itself) is read from CIF as None
        return [str(getattr(self, f"site_symmetry_{i:}", None) or ".")
                for i in range(1, 3 + 1)]


class RestrAngleL(LoopN):
    """Angle restraints of the crystal."""

    ITEM_CLASS = RestrAngle
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(RestrAngleL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class RestrTorsion(ItemN):
    """Restraint of the torsion angle (degrees, IUPAC sign) of four atoms.

``target`` and ``target_weight_param`` (its standard uncertainty sigma)
are given by the user and never recalculated; the refinement adds
((value - target) / sigma)^2 to chi_sq. ``diff`` = value - target is
calculated (Crystal.update_restraints, after every refinement). The atoms
are equivalent positions of atom_site atoms: the symmetry operation n of
the site symmetry code 'n_klm' (as in geom_bond) is kept, the lattice
translation is the one that puts the atom closest to the previous atom of
the chain (to the middle atom for an angle).
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "atom_site_label_4", "target", "target_weight_param",)
    ATTR_MANDATORY_TYPES = (str, str, str, str,  float, float,)
    ATTR_MANDATORY_CIF = ("atom_site_label_1", "atom_site_label_2", "atom_site_label_3", "atom_site_label_4", "target", "target_weight_param",)

    ATTR_OPTIONAL_NAMES = ("site_symmetry_1", "site_symmetry_2", "site_symmetry_3", "site_symmetry_4", "diff", "details",)
    ATTR_OPTIONAL_TYPES = (str, str, str, str,  float, str,)
    ATTR_OPTIONAL_CIF = ("site_symmetry_1", "site_symmetry_2", "site_symmetry_3", "site_symmetry_4", "diff", "details",)

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

    D_FORMATS = {"target": "{:.2f}", "target_weight_param": "{:.2f}", "diff": "{:.2f}"}
    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "restr_torsion"

    def __init__(self, **kwargs) -> NoReturn:
        super(RestrTorsion, self).__init__()
        self.__dict__["D_MIN"] = {"target_weight_param": 0.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_labels(self) -> list:
        return [getattr(self, f"atom_site_label_{i:}") for i in range(1, 4 + 1)]

    def get_site_symmetry(self) -> list:
        # '.' (the atom itself) is read from CIF as None
        return [str(getattr(self, f"site_symmetry_{i:}", None) or ".")
                for i in range(1, 4 + 1)]


class RestrTorsionL(LoopN):
    """Torsion angle restraints of the crystal."""

    ITEM_CLASS = RestrTorsion
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(RestrTorsionL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
