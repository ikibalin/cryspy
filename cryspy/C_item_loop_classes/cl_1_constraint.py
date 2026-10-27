"""Classes Constraint, ConstraintL: hard and soft constraints between parameters."""
from typing import Dict, List, NoReturn, Tuple

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.linear_constraints import \
    parse_linear_constraint, parse_soft_constraint

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Constraint(ItemN):
    """Constraint between refined parameters, hard (x+y=1) or soft.

    The parameters are named by their marks: the letters written after
    the parentheses of a refined value in rcif (``0.7()x``). Without sigma
    the expression is a linear equation of the marks that holds exactly,
    e.g. ``x+y=1`` keeps the sum of two occupancies equal to 1 during the
    refinement. With sigma the constraint is soft (a restraint): any
    equation ``f = g`` of the marks adds ``((f-g)/sigma)**2`` to chi-square.

    Attributes
    ----------
        - expression (mandatory)
        - id, sigma (optional)
    """

    ATTR_MANDATORY_NAMES = ("expression", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("expression", )

    ATTR_OPTIONAL_NAMES = ("id", "sigma")
    ATTR_OPTIONAL_TYPES = (str, float)
    ATTR_OPTIONAL_CIF = ("id", "sigma")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}

    PREFIX = "constraint"

    def __init__(self, **kwargs) -> NoReturn:
        super(Constraint, self).__init__()

        self.__dict__["D_MIN"] = {"sigma": 0.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def is_soft(self) -> bool:
        """True for a restraint: sigma is given."""
        return self.is_attribute("sigma") and self.sigma is not None

    def get_linear_equation(self) -> Tuple[Dict[str, float], float]:
        """Coefficients of the marks and the right-hand side of the equation."""
        return parse_linear_constraint(self.expression)


class ConstraintL(LoopN):
    """Constraints between refined parameters (x+y=1).

    Each expression is written with the marks of the refined values
    (``0.7()x``). Without sigma the constraint holds exactly (linear
    equations only); with sigma it is soft. Example::

        loop_
        _constraint_id
        _constraint_expression
        _constraint_sigma
        1 "x+y=1"       .
        2 "a-2*b=0"     .
        3 "sqrt(a**2+b**2)=2.5"  0.05
    """

    ITEM_CLASS = Constraint
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ConstraintL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_linear_equations(self) -> List[Tuple[Dict[str, float], float]]:
        """Coefficients of the marks and the right-hand side of each hard
        constraint."""
        return [item.get_linear_equation() for item in self.items if not item.is_soft()]

    def get_soft_constraints(self) -> List[Tuple[str, float, set]]:
        """Expression, sigma and marks of each soft constraint."""
        return [(item.expression, float(item.sigma), parse_soft_constraint(item.expression)[1])
                for item in self.items if item.is_soft()]
