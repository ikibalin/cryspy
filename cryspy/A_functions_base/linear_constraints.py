"""Linear constraints between refined parameters (``x + y = 2``).

A constraint is written with the marks of the refined values, the letters
after the parentheses in rcif (``0.7()x``). It is a linear equation:
``sum_i c_i * p_i = b``.
"""
import ast
from typing import Dict, List, Tuple

import numpy


def _parse_linear_node(node) -> Tuple[Dict[str, float], float]:
    """Linear form of an expression node: (coefficients of marks, constant)."""
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) \
            and not isinstance(node.value, bool):
        return {}, float(node.value)
    if isinstance(node, ast.Name):
        return {node.id: 1.}, 0.
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        d_coeff, const = _parse_linear_node(node.operand)
        if isinstance(node.op, ast.USub):
            return {key: -val for key, val in d_coeff.items()}, -const
        return d_coeff, const
    if isinstance(node, ast.BinOp):
        d_left, c_left = _parse_linear_node(node.left)
        d_right, c_right = _parse_linear_node(node.right)
        if isinstance(node.op, (ast.Add, ast.Sub)):
            sign = 1. if isinstance(node.op, ast.Add) else -1.
            d_coeff = dict(d_left)
            for key, val in d_right.items():
                d_coeff[key] = d_coeff.get(key, 0.) + sign * val
            return d_coeff, c_left + sign * c_right
        if isinstance(node.op, ast.Mult):
            if len(d_left) == 0:
                return {key: c_left * val for key, val in d_right.items()}, c_left * c_right
            if len(d_right) == 0:
                return {key: c_right * val for key, val in d_left.items()}, c_left * c_right
            raise ValueError("a product of two parameters is not linear")
        if isinstance(node.op, ast.Div):
            if len(d_right) != 0:
                raise ValueError("division by a parameter is not linear")
            if c_right == 0.:
                raise ValueError("division by zero")
            return {key: val / c_right for key, val in d_left.items()}, c_left / c_right
    raise ValueError(f"'{ast.unparse(node):}' is not allowed in a linear constraint")


def parse_linear_constraint(expression: str) -> Tuple[Dict[str, float], float]:
    """Parse a linear constraint written with the marks of parameters.

    Parameters
    ----------
    expression : str
        Equation such as ``"x + y = 2"``, ``"2*a - b = 0"`` or
        ``"occ_fe = 1 - occ_co"``. Numbers, marks, ``+``, ``-``, ``*`` and
        ``/`` by numbers and parentheses are allowed. Without ``=`` the
        expression is equal to zero.

    Returns
    -------
    d_coeff : dict
        Coefficient of each mark (``{"x": 1., "y": 1.}``); marks whose
        coefficients cancel are left out.
    const : float
        Right-hand side ``b`` of ``sum_i c_i * p_i = b``.

    Raises
    ------
    ValueError
        The expression is not a linear equation of the marks.
    """
    l_part = expression.split("=")
    if len(l_part) > 2:
        raise ValueError(f"constraint '{expression:}' has more than one '='")
    d_coeff, const = {}, 0.
    for sign, part in zip((1., -1.), l_part):
        if part.strip() == "":
            raise ValueError(f"constraint '{expression:}' has an empty side")
        try:
            node = ast.parse(part.strip(), mode="eval").body
        except SyntaxError:
            raise ValueError(f"constraint '{expression:}' cannot be read") from None
        d_part, c_part = _parse_linear_node(node)
        for key, val in d_part.items():
            d_coeff[key] = d_coeff.get(key, 0.) + sign * val
        const -= sign * c_part
    d_coeff = {key: val for key, val in d_coeff.items() if val != 0.}
    if len(d_coeff) == 0:
        raise ValueError(f"constraint '{expression:}' has no parameters")
    return d_coeff, const


def calc_dependent_parameters(matrix_a, vector_b, tolerance: float = 1e-10):
    """Split the parameters of the linear system ``A p = b`` into free and
    dependent ones: ``p_dep = vector_c + matrix_m p_free``.

    The last parameters in the order of the columns are chosen as the
    dependent ones, so the parameters written first in the file stay free.

    Parameters
    ----------
    matrix_a : numpy.ndarray
        Coefficients, shape (n_equation, n_parameter).
    vector_b : numpy.ndarray
        Right-hand sides, shape (n_equation,).
    tolerance : float
        Values smaller than it are zero.

    Returns
    -------
    index_free, index_dep : list of int
        Columns of the free and of the dependent parameters.
    matrix_m : numpy.ndarray
        Shape (n_dep, n_free).
    vector_c : numpy.ndarray
        Shape (n_dep,).

    Raises
    ------
    ValueError
        The equations contradict each other.
    """
    matrix = numpy.array(matrix_a, dtype=float)[:, ::-1]
    vector = numpy.array(vector_b, dtype=float).copy()
    n_eq, n_p = matrix.shape
    l_pivot = []
    i_row = 0
    for i_col in range(n_p):
        if i_row >= n_eq:
            break
        i_max = i_row + int(numpy.argmax(numpy.abs(matrix[i_row:, i_col])))
        if abs(matrix[i_max, i_col]) <= tolerance:
            continue
        matrix[[i_row, i_max]] = matrix[[i_max, i_row]]
        vector[[i_row, i_max]] = vector[[i_max, i_row]]
        coeff = matrix[i_row, i_col]
        matrix[i_row] /= coeff
        vector[i_row] /= coeff
        for i_other in range(n_eq):
            if i_other != i_row:
                factor = matrix[i_other, i_col]
                matrix[i_other] -= factor * matrix[i_row]
                vector[i_other] -= factor * vector[i_row]
        l_pivot.append(i_col)
        i_row += 1
    if numpy.any(numpy.abs(vector[i_row:]) > tolerance * max(1., numpy.abs(vector_b).max(initial=0.))):
        raise ValueError("the constraints contradict each other")
    l_free_rev = [i_col for i_col in range(n_p) if i_col not in l_pivot]
    index_dep = [n_p - 1 - i_col for i_col in l_pivot]
    index_free = sorted(n_p - 1 - i_col for i_col in l_free_rev)
    l_free_rev_sorted = [n_p - 1 - i_free for i_free in index_free]
    matrix_m = -matrix[:len(l_pivot)][:, l_free_rev_sorted]
    vector_c = vector[:len(l_pivot)]
    return index_free, index_dep, matrix_m, vector_c


def form_linear_system(
        l_equation: List[Tuple[List[Tuple[float, object]], float]],
        parameter_names: List[object]):
    """Matrix and right-hand side of ``A p = b`` for equations given as
    ``([(coeff, parameter_name), ...], b)``; parameter_names give the
    order of the columns."""
    matrix_a = numpy.zeros((len(l_equation), len(parameter_names)), dtype=float)
    vector_b = numpy.zeros((len(l_equation),), dtype=float)
    for i_eq, (l_term, const) in enumerate(l_equation):
        for coeff, name in l_term:
            matrix_a[i_eq, parameter_names.index(name)] += coeff
        vector_b[i_eq] = const
    return matrix_a, vector_b
