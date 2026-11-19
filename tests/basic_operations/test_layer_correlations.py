"""Layers without correlations between them: Warren line shape and layer search."""
import math

import numpy
import pytest
from scipy.integrate import quad

import cryspy
from cryspy.A_functions_base.layer_correlations import (
    calc_layer_cuts,
    calc_layers_of_atoms,
    calc_rod_powder_factor,
    calc_rods,
    calc_warren_function,
)
from cryspy.procedure_rhochi.rhochi_layer_correlation import parse_normal

TBTE3 = """
data_tbte3
_cell_length_a 4.2802
_cell_length_b 25.2855
_cell_length_c 4.2903
_cell_angle_alpha 90.
_cell_angle_beta 90.
_cell_angle_gamma 90.
_space_group_name_H-M_alt "C m c m"
_space_group_IT_coordinate_system_code abc
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
Tb  Tb 0. 0.170055 0.25 1. Biso 0.16
Te1 Te 0. 0.294884 0.25 1. Biso 0.37
Te2 Te 0. 0.571669 0.25 1. Biso 0.79
Te3 Te 0. 0.928438 0.25 1. Biso 0.51
"""

ORTHO = numpy.array([4., 4., 6., 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])


def warren_by_quad(a: float) -> float:
    x_peak = math.sqrt(max(a, 0.))
    def func(x):
        return math.exp(-(x * x - a) ** 2)
    return quad(func, 0., x_peak, limit=500)[0] + quad(func, x_peak, x_peak + 4., limit=500)[0]


def test_warren_function_special_values():
    # F(0) = int exp(-x^4) dx = Gamma(5/4)
    assert calc_warren_function(numpy.array([0.]))[0] == pytest.approx(math.gamma(1.25), rel=1e-10)
    # the table (linear interpolation) and the Gauss-Hermite rule agree
    # with a direct integration
    a = numpy.array([-4., -1., 0.3, 2., 9.7, 29.9, 30.1, 120.])
    expected = numpy.array([warren_by_quad(value) for value in a])
    assert numpy.allclose(calc_warren_function(a), expected, rtol=1e-6, atol=1e-14)
    # the tail of the line: sqrt(pi) / (2 sqrt(a))
    assert calc_warren_function(numpy.array([1e4]))[0] == pytest.approx(
        math.sqrt(math.pi) / 200., rel=1e-6)


@pytest.mark.parametrize("s", [0.395, 0.4, 0.404, 0.42, 0.5, 0.7])
def test_rod_powder_factor_against_integral(s):
    # J = int_0^s p(v - g) / sqrt(s^2 - v^2) dv, p(x) = L exp(-pi L^2 x^2):
    # the formula with sqrt(s + v) ~ sqrt(s + g) is within 1 % at the front
    # and exact in the tail
    g, length = 0.4, 80.
    expected = quad(lambda v: length * math.exp(-math.pi * (length * (v - g)) ** 2)
                    / math.sqrt(s * s - v * v), 0., s, limit=800,
                    points=[g] if g < s else None)[0]
    value = calc_rod_powder_factor(numpy.array([s]), numpy.array([g]), length)[0]
    assert value == pytest.approx(expected, rel=1e-2)


def test_layer_cuts_and_layers():
    # two equivalent layers at y = 0 and 1/2 of a 12 A long cell
    heights = numpy.array([0.02, -0.02, 0.52, 0.48])
    cuts, gaps = calc_layer_cuts(heights, 12.)
    assert numpy.allclose(cuts, [0.25, 0.75])
    assert gaps[0] == pytest.approx(0.46 * 12.)
    fract = numpy.array([[0., 0., 0.5, 0.5], [0.02, 0.98, 0.52, 0.48], [0., 0., 0., 0.]])
    fract_layer, layer_index, n_layer = calc_layers_of_atoms(fract, (0, 1, 0), ORTHO)
    assert n_layer == 2
    # the atom at y = 0.98 is moved next to the one at 0.02
    assert numpy.allclose(fract_layer[1], [0.02, -0.02, 0.52, 0.48])
    assert layer_index[0] == layer_index[1] != layer_index[2] == layer_index[3]


def test_rods_of_square_layers():
    index_rod, g_rod, g_ccs = calc_rods((0, 0, 1), ORTHO, 0.36)
    # (10), (01), (-10), (0-1) at 1/4 and (11)... at sqrt(2)/4; no rod at the origin
    assert g_rod.size == 8
    assert numpy.allclose(numpy.sort(g_rod), [0.25] * 4 + [math.sqrt(2.) / 4.] * 4)
    assert numpy.allclose(g_ccs[2], 0.)


def test_parse_normal():
    assert parse_normal("0 2 0") == (0, 1, 0)
    assert parse_normal("(1,0,-1)") == (1, 0, -1)
    assert parse_normal("") is None
    with pytest.raises(ValueError):
        parse_normal("1 0")


def test_suggest_layer_direction_tbte3():
    # TbTe3 (Cmcm): van der Waals gap between Te sheets normal to b
    rcif_object = cryspy.str_to_globaln(TBTE3)
    crystal = rcif_object.crystal_tbte3
    from cryspy.procedure_rhochi.rhochi_layer_correlation import calc_layer_candidates
    best = calc_layer_candidates(crystal, magnetic_labels=["Tb"])[0]
    assert best["normal"] == (0, 1, 0)
    assert best["n_layer"] == 2
    assert best["ratio"] > 1.3
    assert best["m_inter"] / best["m_intra"] > 2.
    text = cryspy.suggest_layer_direction(rcif_object)
    assert "planes (0 1 0)" in text


def test_assistant_finds_procedures_and_recipe():
    from cryspy.procedure_assistant.methods import form_matching_methods_text
    from cryspy.procedure_assistant.procedures import describe_procedure
    from cryspy.procedure_assistant.recipes import take_recipe
    rcif_object = cryspy.str_to_globaln(TBTE3)
    text = form_matching_methods_text(rcif_object, "Add 2D correlations")
    assert "- add_layer_correlation:" in text
    text = form_matching_methods_text(rcif_object, "Suggest the direction of the layers")
    assert "- suggest_layer_direction:" in text
    for function in (cryspy.suggest_layer_direction, cryspy.add_layer_correlation,
                     cryspy.report_layer_correlation):
        description = describe_procedure(function)
        assert all(parameter.get("question", "x") for parameter in description["parameters"][1:])
    assert take_recipe("layer_correlation") is not None
    assert "Warren" in cryspy.get_help("sawtooth asymmetric peak 2D correlations")
