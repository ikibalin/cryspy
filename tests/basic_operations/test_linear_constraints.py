"""Linear constraints between refined parameters (x+y=1)."""
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.linear_constraints import (
    calc_dependent_parameters, calc_soft_constraint_chi_sq, parse_linear_constraint,
    parse_soft_constraint)

F_NAME = os.path.join(os.path.dirname(__file__), "..", "functional", "twin_tbte3", "main.rcif")
PATH_1 = "tbte3.atom_site[Te1].fract_y"
PATH_2 = "tbte3.atom_site[Te2].fract_y"


def test_parse_linear_constraint():
    assert parse_linear_constraint("x + y = 2") == ({"x": 1., "y": 1.}, 2.)
    assert parse_linear_constraint("2*(a - b)/4 = 1 - c") == ({"a": 0.5, "b": -0.5, "c": 1.}, 1.)
    assert parse_linear_constraint("x - y") == ({"x": 1., "y": -1.}, 0.)
    for expression in ("x*y = 1", "1/x = 2", "x = y = z", "2 = 2", "x + = 1", "sin(x) = 0"):
        with pytest.raises(ValueError):
            parse_linear_constraint(expression)


def test_calc_dependent_parameters():
    # p0 + p1 = 1, p1 - 2*p2 = 0: the last parameters depend on p0
    matrix_a = numpy.array([[1., 1., 0.], [0., 1., -2.]])
    vector_b = numpy.array([1., 0.])
    index_free, index_dep, matrix_m, vector_c = calc_dependent_parameters(matrix_a, vector_b)
    assert index_free == [0]
    assert sorted(index_dep) == [1, 2]
    p_free = numpy.array([0.3])
    p_full = numpy.zeros(3)
    p_full[index_free] = p_free
    p_full[index_dep] = vector_c + matrix_m.dot(p_free)
    assert numpy.allclose(matrix_a.dot(p_full), vector_b)
    with pytest.raises(ValueError):
        calc_dependent_parameters(numpy.array([[1., 1.], [1., 1.]]), numpy.array([1., 2.]))


def test_add_constraint_and_dictionary():
    obj = cryspy.load_file(F_NAME)
    cryspy.add_constraint(obj, "y1 + y2 = 1.5", {"y1": PATH_1, "y2": PATH_2})
    assert cryspy.take_constraints(obj) == ["y1 + y2 = 1.5"]
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    assert isinstance(obj_2.items[0], cryspy.ConstraintL)
    l_term, const = obj_2.get_dictionary()["linear_equations"][0]
    assert const == 1.5
    assert [coeff for coeff, name in l_term] == [1., 1.]
    d_marked = {d_par["mark"]: d_par["path"] for d_par in cryspy.take_marked_parameters(obj_2)
                if d_par["mark"] != ""}
    assert d_marked == {"y1": "crystal_tbte3.atom_site[Te1].fract_y",
                        "y2": "crystal_tbte3.atom_site[Te2].fract_y"}
    assert cryspy.propose_mark(obj_2, PATH_1) == "y3"


def test_add_constraint_errors():
    obj = cryspy.load_file(F_NAME)
    with pytest.raises(cryspy.RcifFormatError):
        cryspy.add_constraint(obj, "z = 1")
    with pytest.raises(cryspy.RcifFormatError):
        cryspy.add_constraint(obj, "a*b = 1", {"a": PATH_1, "b": PATH_2})
    cryspy.add_constraint(obj, "a + b = 1.5", {"a": PATH_1, "b": PATH_2})
    with pytest.raises(cryspy.RcifFormatError):
        # the mark is used by another parameter
        cryspy.set_parameter_marks(obj, {"a": PATH_2})
    assert cryspy.remove_constraint(obj, "1")
    assert cryspy.take_constraints(obj) == []


def test_unknown_mark_in_constraint_loop():
    s_cont = open(F_NAME).read().replace("global_", 'global_\n\nloop_\n_constraint_id\n'
                                         '_constraint_expression\n1 "q = 1"\n', 1)
    obj = cryspy.str_to_globaln(s_cont)
    with pytest.raises(cryspy.RcifFormatError):
        obj.get_dictionary()


def test_assistant_tools():
    obj = cryspy.load_file(F_NAME)
    text = cryspy.call_assistant_tool(
        obj, "add_constraint",
        {"expression": "y1+y2=1.5", "parameters": f"y1={PATH_1}, y2={PATH_2}"})["text"]
    assert "Constraint 1 added" in text
    text = cryspy.call_assistant_tool(obj, "list_constraints", {})["text"]
    assert "y1+y2=1.5" in text and "y2 = crystal_tbte3.atom_site[Te2].fract_y" in text


def test_soft_constraint_expression():
    node, marks = parse_soft_constraint("sqrt(a**2 + b**2) = 2.5")
    assert marks == {"a", "b"}
    assert calc_soft_constraint_chi_sq("sqrt(a**2 + b**2) = 2.5", 0.1, {"a": 1.5, "b": 2.1}) == \
        pytest.approx(((numpy.hypot(1.5, 2.1) - 2.5) / 0.1) ** 2)
    assert calc_soft_constraint_chi_sq("cosd(t) = 0.5", 0.01, {"t": 60.}) == pytest.approx(0., abs=1e-20)
    for expression in ("__import__('os')", "a.real = 1", "sqrt = 1", "max(a, b) = 1", "1 = 2"):
        with pytest.raises(ValueError):
            parse_soft_constraint(expression)


def test_add_soft_constraint_and_dictionary():
    obj = cryspy.load_file(F_NAME)
    cryspy.add_constraint(obj, "y1 + y2 = 1.5", {"y1": PATH_1, "y2": PATH_2}, sigma=0.01)
    cryspy.add_constraint(obj, "y1 - y2 = 0.36")
    assert cryspy.take_constraints(obj) == ["y1 + y2 = 1.5 (sigma 0.01)", "y1 - y2 = 0.36"]
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    d_global = obj_2.get_dictionary()
    assert len(d_global["linear_equations"]) == 1
    expression, sigma, d_mark_way = d_global["soft_constraints"][0]
    assert (expression, sigma, sorted(d_mark_way)) == ("y1 + y2 = 1.5", 0.01, ["y1", "y2"])
    with pytest.raises(cryspy.RcifFormatError):
        cryspy.add_constraint(obj, "y1 = 1", sigma=0.)
    text = cryspy.call_assistant_tool(
        obj, "add_constraint", {"expression": "y1**2 = 0.86", "sigma": 0.001})["text"]
    assert "(sigma 0.001)" in text
