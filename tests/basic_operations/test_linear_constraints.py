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


def test_marked_parameters_are_found_in_dictionary():
    """Occupancies and cell lengths are not in the table of
    transfer_cif_names_to_dict_names: their ways are found in the
    dictionary of the data block (the refinement failed with
    'The constrained parameter ... is not refined')."""
    s_cont = open(F_NAME).read()
    s_cont = s_cont.replace("Te1     Te   0.000000   0.928930   0.250000   1.0",
                            "Te1     Te   0.000000   0.928930   0.250000   0.79()occ1")
    s_cont = s_cont.replace("Te2     Te   0.000000   0.571030   0.250000   1.0",
                            "Te2     Te   0.000000   0.571030   0.250000   1.0()occ2")
    s_cont = s_cont.replace("_cell_length_c 4.290300", "_cell_length_c 4.290300()c")
    s_cont = s_cont.replace("global_", 'global_\n\nloop_\n_constraint_id\n_constraint_expression\n'
                            '1 "occ1+occ2=1"\n2 "c=4.3"\n', 1)
    obj = cryspy.str_to_globaln(s_cont)
    d_global = obj.get_dictionary()
    assert d_global["linear_equations"] == [
        ([(1., ("crystal_tbte3", "atom_occupancy", (1, ))),
          (1., ("crystal_tbte3", "atom_occupancy", (2, )))], 1.),
        ([(1., ("crystal_tbte3", "unit_cell_parameters", (2, )))], 4.3)]
    # the object is not changed by the search
    assert obj.crystal_tbte3.atom_site["Te2"].occupancy == 1.
    assert obj.crystal_tbte3.cell.length_c == 4.2903


U_1 = "tbte3.atom_site[Te1].u_iso_or_equiv"
U_2 = "tbte3.atom_site[Te2].u_iso_or_equiv"
COEFF_B = 8. * numpy.pi ** 2


def take_u_constraint_object(s_replace: tuple = ()):
    s_cont = open(F_NAME).read()
    for old, new in s_replace:
        s_cont = s_cont.replace(old, new)
    obj = cryspy.str_to_globaln(s_cont)
    obj.fix_variables()
    cryspy.add_constraint(obj, "u1 + u2 = 0.02", {"u1": U_1, "u2": U_2})
    return obj


def test_constraint_in_units_of_the_object():
    """The dictionary keeps b_iso = 8 pi^2 u_iso: the constraint
    u1 + u2 = 0.02 is b1 + b2 = 8 pi^2 0.02, not b1 + b2 = 0.02."""
    obj = take_u_constraint_object()
    l_term, const = obj.get_dictionary()["linear_equations"][0]
    assert [name for coeff, name in l_term] == [
        ("crystal_tbte3", "atom_b_iso", (1, )), ("crystal_tbte3", "atom_b_iso", (2, ))]
    b_1, b_2 = 0.012, 0.008  # u1 + u2 = 0.02 in b
    assert sum(coeff * b for (coeff, name), b in zip(l_term, (b_1 * COEFF_B, b_2 * COEFF_B))) == \
        pytest.approx(const)
    assert const == pytest.approx(0.02)


def test_cell_angle_constraint_in_degrees():
    """Cell angles are in radians in the dictionary."""
    s_cont = open(F_NAME).read().replace("C m c m", "P -1").replace(
        "_cell_angle_alpha 90.000000", "_cell_angle_alpha 90.000000()al").replace(
        "_cell_angle_beta 90.000000", "_cell_angle_beta 90.000000()be")
    s_cont = s_cont.replace("global_", 'global_\n\nloop_\n_constraint_id\n_constraint_expression\n'
                            '1 "al+be=181"\n', 1)
    obj = cryspy.str_to_globaln(s_cont)
    l_term, const = obj.get_dictionary()["linear_equations"][0]
    assert [name for coeff, name in l_term] == [
        ("crystal_tbte3", "unit_cell_parameters", (3, )),
        ("crystal_tbte3", "unit_cell_parameters", (4, ))]
    assert [coeff for coeff, name in l_term] == pytest.approx([180. / numpy.pi] * 2)
    assert const == pytest.approx(181.)


def test_not_refined_constrained_parameter_is_refined():
    """x+y=1 with only x refined: y is refined too (the refinement failed
    with 'The constrained parameter ... is not refined')."""
    obj = take_u_constraint_object()
    atom_site = obj.crystal_tbte3.atom_site
    atom_site["Te2"].u_iso_or_equiv_refinement = False
    d_global = obj.get_dictionary()
    assert atom_site["Te2"].u_iso_or_equiv_refinement
    assert len(d_global["linear_equations"]) == 1
    # nothing of the constraint is refined: it is not used
    atom_site["Te1"].u_iso_or_equiv_refinement = False
    atom_site["Te2"].u_iso_or_equiv_refinement = False
    assert "linear_equations" not in obj.get_dictionary()
    assert not atom_site["Te1"].u_iso_or_equiv_refinement


def test_equal_marks_refine_both_parameters():
    s_cont = open(F_NAME).read()
    s_cont = s_cont.replace("Te1     Te   0.000000   0.928930   0.250000   1.0   Uiso   0.0",
                            "Te1     Te   0.000000   0.928930   0.250000   1.0   Uiso   0.01()a")
    s_cont = s_cont.replace("Te2     Te   0.000000   0.571030   0.250000   1.0   Uiso   0.0",
                            "Te2     Te   0.000000   0.571030   0.250000   1.0   Uiso   0.02()a")
    obj = cryspy.str_to_globaln(s_cont)
    obj.crystal_tbte3.atom_site["Te2"].u_iso_or_equiv_refinement = False
    d_global = obj.get_dictionary()
    assert obj.crystal_tbte3.atom_site["Te2"].u_iso_or_equiv_refinement
    assert d_global["linear_constraints"] == [
        [(pytest.approx(1. / COEFF_B), ("crystal_tbte3", "atom_b_iso", (1, ))),
         (pytest.approx(-1. / COEFF_B), ("crystal_tbte3", "atom_b_iso", (2, )))]]


def test_soft_constraint_in_units_of_the_object():
    from cryspy.procedure_rhochi.rhochi_by_dictionary import calc_punishement_function
    obj = take_u_constraint_object()
    cryspy.remove_constraint(obj, "1")
    cryspy.add_constraint(obj, "u1 + u2 = 0.02", sigma=0.001)
    atom_site = obj.crystal_tbte3.atom_site
    atom_site["Te1"].u_iso_or_equiv = 0.011
    atom_site["Te2"].u_iso_or_equiv = 0.010
    assert calc_punishement_function(obj.get_dictionary()) == pytest.approx(1.)


def test_u_iso_without_adp_type():
    """Without _atom_site_adp_type a given U_iso_or_equiv is used (Uiso)."""
    obj = take_u_constraint_object((("_atom_site_adp_type\n", ""), ("   Uiso   ", "   ")))
    atom_site = obj.crystal_tbte3.atom_site
    atom_site["Te1"].u_iso_or_equiv = 0.01
    d_global = obj.get_dictionary()
    assert d_global["crystal_tbte3"]["atom_b_iso"][1] == pytest.approx(0.01 * COEFF_B)
    assert d_global["linear_equations"][0][0][0][1] == ("crystal_tbte3", "atom_b_iso", (1, ))


def test_constrained_parameter_not_used_by_the_model():
    """U of an atom with adp_type Biso does not change the calculation."""
    obj = take_u_constraint_object((("   Uiso   ", "   Biso   "), ))
    for atom in obj.crystal_tbte3.atom_site.items:
        atom.b_iso_or_equiv = 0.5
    with pytest.raises(cryspy.RcifFormatError, match="adp_type"):
        obj.get_dictionary()


def test_marks_in_data_block_with_capital_letters():
    """data_Pbcn: names of get_variable_names are in lower case, the
    parameters are found and changed by them (the calculation failed with
    'cannot unpack non-iterable NoneType object')."""
    s_cont = open(F_NAME).read().replace("data_tbte3", "data_TbTe3")
    s_cont = s_cont.replace("Te1     Te   0.000000   0.928930   0.250000   1.0   Uiso   0.0",
                            "Te1     Te   0.000000   0.928930   0.250000   1.0   Uiso   0.01()rr")
    s_cont = s_cont.replace("Te2     Te   0.000000   0.571030   0.250000   1.0   Uiso   0.0",
                            "Te2     Te   0.000000   0.571030   0.250000   1.0   Uiso   0.01()rr")
    obj = cryspy.str_to_globaln(s_cont)
    d_global = obj.get_dictionary()
    assert d_global["linear_constraints"] == [
        [(pytest.approx(1. / COEFF_B), ("crystal_tbte3", "atom_b_iso", (1, ))),
         (pytest.approx(-1. / COEFF_B), ("crystal_tbte3", "atom_b_iso", (2, )))]]
    name = (("global", ""), ("crystal", "tbte3"), ("atom_site", ""), ("u_iso_or_equiv", 1))
    obj.set_variable_by_name(name, 0.02)
    assert obj.get_variable_by_name(name) == 0.02
