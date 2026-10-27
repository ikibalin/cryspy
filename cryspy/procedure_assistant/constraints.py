"""Linear constraints between refined parameters (x+y=1) of a CrysPy object.

A constraint names the parameters by their marks, the letters written after
the parentheses of a refined value in rcif (``0.7()x``). The constraints are
kept in the loop ``constraint`` of the global object:

    loop_
    _constraint_id
    _constraint_expression
    1 "x+y=1"

CrysPy Editor, AI assistants and scripts create constraints through
add_constraint, so the marks and the loop are always consistent.
"""
import re
from typing import Dict, List, Optional

from cryspy.A_functions_base.linear_constraints import (
    parse_linear_constraint, parse_soft_constraint)
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_constraint import Constraint, ConstraintL
from cryspy.exceptions import RcifFormatError

from cryspy.procedure_assistant.parameters import (
    ParameterRecord, take_parameter_record, take_parameter_records)

RE_MARK = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def take_mark(record: ParameterRecord) -> str:
    """Mark of the parameter ('' when it has none)."""
    return str(getattr(record.item, f"{record.attribute:}_mark", ""))


def set_mark(record: ParameterRecord, mark: str):
    """Write the mark of the parameter."""
    setattr(record.item, f"{record.attribute:}_mark", mark)


def take_constraint_loop(cryspy_object: GlobalN) -> Optional[ConstraintL]:
    for item in cryspy_object.items:
        if isinstance(item, ConstraintL):
            return item
    return None


def take_constraints(cryspy_object) -> List[str]:
    """Expressions of the constraints of the object; a soft one ends with
    its sigma: 'x+y=1 (sigma 0.01)'."""
    if not isinstance(cryspy_object, GlobalN):
        return []
    loop = take_constraint_loop(cryspy_object)
    if loop is None:
        return []
    return [str(item.expression) + (f" (sigma {item.sigma:g})" if item.is_soft() else "")
            for item in loop.items]


def take_marked_parameters(cryspy_object) -> List[dict]:
    """Refinable parameters with their marks, for a dialog or a model.

    Returns
    -------
    list of dict
        Keys "path", "mark", "value", "refined" and "constrained" (fixed
        by symmetry, cannot be used).
    """
    l_parameter = []
    for record in take_parameter_records(cryspy_object):
        l_parameter.append({
            "path": record.path, "mark": take_mark(record), "value": record.value,
            "refined": record.flag_refinement, "constrained": record.flag_constraint})
    return l_parameter


def propose_mark(cryspy_object, path: str = "") -> str:
    """Mark not used in the object yet, made from the attribute name of
    the path ('fract_y' -> 'y1', 'occupancy' -> 'occ1') or 'p1', 'p2', ..."""
    l_used = {take_mark(record).rstrip("-+") for record in take_parameter_records(cryspy_object)}
    l_word = [word for word in re.split(r"[^A-Za-z0-9]+", path.split(".")[-1]) if word != ""]
    if len(l_word) == 0:
        base = ""
    elif len(l_word[-1]) <= 2:
        base = l_word[-1].lower()
    else:
        base = l_word[0][:3].lower()
    if base == "" or not base[0].isalpha():
        base = "p"
    index = 1
    while f"{base:}{index:}" in l_used:
        index += 1
    return f"{base:}{index:}"


def check_mark(mark: str) -> str:
    mark_sh = mark.strip()
    if RE_MARK.match(mark_sh) is None:
        raise RcifFormatError(
            f"'{mark:}' cannot be a mark: use a letter followed by letters, "
            "digits or '_' (x, occ1, a_2)")
    return mark_sh


def set_parameter_marks(cryspy_object, d_mark_path: Dict[str, str]) -> List[ParameterRecord]:
    """Give marks to parameters and switch on their refinement.

    Parameters
    ----------
    d_mark_path : dict
        Mark -> path of the parameter ('x' -> 'tbte3.atom_site[Te1].fract_y').

    Raises
    ------
    RcifFormatError
        A path is unknown, a parameter is fixed by symmetry or a mark is
        not a name.
    """
    l_record = []
    for mark, path in d_mark_path.items():
        mark_sh = check_mark(mark)
        record = take_parameter_record(cryspy_object, path)
        if record is None:
            raise RcifFormatError(f"There is no refinable parameter '{path:}'.")
        if record.flag_constraint:
            raise RcifFormatError(
                f"{record.path:} is fixed by symmetry and cannot be constrained.")
        for other in take_parameter_records(cryspy_object):
            if other.path != record.path and take_mark(other).rstrip("-+") == mark_sh:
                raise RcifFormatError(
                    f"The mark {mark_sh:} is used by {other.path:} already: "
                    "choose another mark (equal marks make values equal).")
        l_record.append((mark_sh, record))
    for mark_sh, record in l_record:
        set_mark(record, mark_sh)
        record.set_refinement(True)
    return [record for mark_sh, record in l_record]


def add_constraint(cryspy_object: GlobalN, expression: str,
                   d_mark_path: Optional[Dict[str, str]] = None,
                   sigma: Optional[float] = None) -> Constraint:
    """Add a constraint 'x+y=1' to the loop of constraints.

    Without sigma the constraint is hard: a linear equation that holds
    exactly. With sigma it is soft (a restraint): any equation f = g of
    the marks (``sqrt(a**2+b**2)=2.5``) adds ((f-g)/sigma)**2 to chi-square.

    Parameters
    ----------
    cryspy_object : GlobalN
        Object with the crystals and the experiments.
    expression : str
        Equation of the marks of refined parameters.
    d_mark_path : dict, optional
        Marks to give first: mark -> path of the parameter.
    sigma : float, optional
        Allowed deviation of the soft constraint, in the units of the
        expression; None for a hard constraint.

    Returns
    -------
    Constraint
        The new item of the loop of constraints.

    Raises
    ------
    RcifFormatError
        The expression is not linear or uses a mark of no refined value.
    """
    if not isinstance(cryspy_object, GlobalN):
        raise RcifFormatError("Constraints are kept in the global object only.")
    try:
        if sigma is None:
            l_label = list(parse_linear_constraint(expression)[0])
        else:
            sigma = float(sigma)
            if not sigma > 0.:
                raise ValueError("sigma must be positive")
            l_label = sorted(parse_soft_constraint(expression)[1])
    except ValueError as error:
        raise RcifFormatError(f"Constraint '{expression:}': {error:}.") from None
    if d_mark_path:
        set_parameter_marks(cryspy_object, d_mark_path)
    l_mark = {take_mark(record).rstrip("-+") for record in take_parameter_records(cryspy_object)
              if record.flag_refinement}
    l_unknown = [mark for mark in l_label if mark not in l_mark]
    if len(l_unknown) > 0:
        raise RcifFormatError(
            f"Constraint '{expression:}': no refined parameter is marked "
            f"{', '.join(l_unknown):}.")
    loop = take_constraint_loop(cryspy_object)
    if loop is None:
        loop = ConstraintL()
        cryspy_object.add_items([loop])
    l_id = [str(item.id) for item in loop.items if item.is_attribute("id")]
    index = len(loop.items) + 1
    while str(index) in l_id:
        index += 1
    item = Constraint(id=str(index), expression=expression.strip())
    if sigma is not None:
        item.sigma = sigma
    loop.items.append(item)
    return item


def remove_constraint(cryspy_object: GlobalN, index: str) -> bool:
    """Remove the constraint with this id; True if it was found."""
    loop = take_constraint_loop(cryspy_object) if isinstance(cryspy_object, GlobalN) else None
    if loop is None:
        return False
    for item in loop.items:
        if item.is_attribute("id") and str(item.id) == str(index):
            loop.items.remove(item)
            if len(loop.items) == 0:
                cryspy_object.items.remove(loop)
            return True
    return False
