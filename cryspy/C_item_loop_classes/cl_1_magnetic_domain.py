"""MagneticDomain and MagneticDomainL classes (category magnetic_domain)."""
from typing import NoReturn

import numpy

from cryspy.A_functions_base.function_1_objects import form_items_by_dictionary
from cryspy.A_functions_base.structure_factor_magnetic_k import parse_operation_xyz
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.exceptions import CrysPyError


def parse_magnetic_operation(operation_xyz: str):
    """Rotation (9,), translation (3,) and time reversal (+1 or -1) of a
    magnetic operation written as in mCIF, ``"-y,x,z+1/2,-1"``; without the
    fourth term the time reversal is +1."""
    l_term = [term.strip() for term in str(operation_xyz).split(",")]
    if len(l_term) == 4:
        try:
            time_reversal = float(l_term[3])
        except ValueError:
            time_reversal = 0.
        if abs(abs(time_reversal) - 1.) > 1e-8:
            raise CrysPyError(
                f"The time reversal of the magnetic operation '{operation_xyz:}' "
                "should be +1 or -1.")
    elif len(l_term) == 3:
        time_reversal = 1.
    else:
        raise CrysPyError(
            f"The magnetic operation '{operation_xyz:}' should be 'x,y,z' or 'x,y,z,+1'.")
    rotation, translation = parse_operation_xyz(",".join(l_term[:3]))
    return rotation, translation, time_reversal


class MagneticDomain(ItemN):
    """Magnetic domain of a single crystal experiment.

The domain is the crystal transformed by a magnetic operation of the
paramagnetic group that is not in the magnetic group, written as in mCIF
(operation_xyz "-y,x,z,+1": rotation, translation and time reversal). The
nuclear structure is the same in all domains. Fraction gives the volume
fraction of the domain; the first domain (identity, always present) has
1 - sum of the fractions. The magnetic domains are combined with the
structural twins of jana_twin: every twin has all the magnetic domains.

Fraction is a refined parameter.
    """

    ATTR_MANDATORY_NAMES = ("id", "operation_xyz", "fraction")
    ATTR_MANDATORY_TYPES = (str, str, float)
    ATTR_MANDATORY_CIF = ("id", "operation_xyz", "fraction")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("fraction", )
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

    PREFIX = "magnetic_domain"

    def __init__(self, **kwargs) -> NoReturn:
        super(MagneticDomain, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"fraction": 0., }

        # defined for ani integer and float parameters
        D_MAX = {"fraction": 1., }

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class MagneticDomainL(LoopN):
    """Magnetic domains of a single crystal experiment.

Every row is one more domain (the first one, the identity, is not written):
the magnetic operation operation_xyz ("-y,x,z,+1", as in mCIF) and the
refined volume fraction. The first domain has 1 - sum of the fractions.
    """
    ITEM_CLASS = MagneticDomain
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MagneticDomainL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        """Fractions (with refinement flags), matrices on hkl R^T (9, n) and
        time reversals (n,) of the domains."""
        res = {}
        l_matrix, l_time_reversal = [], []
        for item in self.items:
            rotation, translation, time_reversal = parse_magnetic_operation(item.operation_xyz)
            # the domain gives at H the intensity of the first one at R^T H
            l_matrix.append(rotation.reshape(3, 3).transpose().flatten())
            l_time_reversal.append(time_reversal)
        res["magnetic_domain_fraction"] = numpy.array(self.fraction, dtype=float)
        res["flags_magnetic_domain_fraction"] = numpy.array(self.fraction_refinement, dtype=bool)
        res["magnetic_domain_matrices"] = numpy.array(l_matrix, dtype=float).transpose()
        res["magnetic_domain_time_reversal"] = numpy.array(l_time_reversal, dtype=float)
        return res
