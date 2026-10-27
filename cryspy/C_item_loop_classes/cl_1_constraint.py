"""Classes Constraint, ConstraintL: linear constraints between parameters."""
from typing import Dict, List, NoReturn, Tuple

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.linear_constraints import \
    parse_linear_constraint

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Constraint(ItemN):
    """Linear constraint between refined parameters.

    The parameters are named by their marks: the letters written after
    the parentheses of a refined value in rcif (``0.7()x``). The
    expression is a linear equation of the marks, e.g. ``x+y=1`` keeps
    the sum of two occupancies equal to 1 during the refinement.

    Attributes
    ----------
        - expression (mandatory)
        - id (optional)
    """

    ATTR_MANDATORY_NAMES = ("expression", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("expression", )

    ATTR_OPTIONAL_NAMES = ("id", )
    ATTR_OPTIONAL_TYPES = (str, )
    ATTR_OPTIONAL_CIF = ("id", )

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

        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_linear_equation(self) -> Tuple[Dict[str, float], float]:
        """Coefficients of the marks and the right-hand side of the equation."""
        return parse_linear_constraint(self.expression)


class ConstraintL(LoopN):
    """Linear constraints between refined parameters (x+y=1).

    Each expression is written with the marks of the refined values
    (``0.7()x``). Example::

        loop_
        _constraint_id
        _constraint_expression
        1 "x+y=1"
        2 "a-2*b=0"
    """

    ITEM_CLASS = Constraint
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ConstraintL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_linear_equations(self) -> List[Tuple[Dict[str, float], float]]:
        """Coefficients of the marks and the right-hand side of each equation."""
        return [item.get_linear_equation() for item in self.items]
