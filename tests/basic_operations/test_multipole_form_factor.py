"""Shared multipole form factor of magnetic atoms (magnetic space groups and
superspace groups) and Crystal.add_multipole_form_factor."""
import numpy
import pytest

import cryspy
from cryspy.A_functions_base.magnetic_form_factor import (
    calc_magnetic_form_factor_with_multipoles_by_dictionary)
from cryspy.D_functions_item_loop.function_1_multipoles import (
    parse_multipole_coefficients, take_orbital_of_ion)
from cryspy.exceptions import CrysPyError

RCIF = """global_

data_cu
_cell_length_a 5.4
_cell_length_b 5.4
_cell_length_c 13.1
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0

_space_group_name_H-M_alt 'P 1'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
Cu Cu2+ 0.0 0.0 0.0 1.0 Biso 0.0
La La3+ 0.0 0.0 0.36 1.0 Biso 0.0

loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_lande
Cu Cu2+ 2.0
"""


def test_parse_multipole_coefficients():
    assert parse_multipole_coefficients("P20 0.1, P2-2=-0.05; coeff_P44 1e-2") == {
        "coeff_P20": 0.1, "coeff_P2m2": -0.05, "coeff_P44": 0.01}
    assert parse_multipole_coefficients("") == {}
    with pytest.raises(CrysPyError):
        parse_multipole_coefficients("P23 0.1")


def test_orbital_of_ion_from_knowledge_base():
    assert take_orbital_of_ion("Cu2+") == "3d"
    assert take_orbital_of_ion("Tb3+") == "4f"


def test_spherical_form_factor_is_kept_without_multipoles():
    form_factor = numpy.arange(6.).reshape(3, 2)
    result = calc_magnetic_form_factor_with_multipoles_by_dictionary(
        {}, numpy.zeros((3, 3)), numpy.eye(3).reshape(9, 1), numpy.array(["A", "B"]),
        form_factor)
    assert result.shape == (3, 1, 2)
    numpy.testing.assert_array_equal(result[:, 0, :], form_factor)


def test_add_multipole_form_factor_to_cu():
    obj = cryspy.str_to_globaln(RCIF)
    crystal = obj.crystal_cu
    report = crystal.add_multipole_form_factor("Cu", coefficients="P22 0.2")
    assert "Ion Cu2+, orbital 3d" in report
    item = crystal.atom_rho_multipole["Cu"]
    assert item.coeff_P00 == pytest.approx(1.) and item.coeff_P22 == pytest.approx(0.2)
    assert list(crystal.atom_local_axes.atom_label) == ["Cu"]
    assert {"DCu_c", "DCu_a"} <= set(crystal.atom_site.label)
    # local axes are the Cartesian axes of CrysPy: x || a*, z || c
    matrix = crystal.get_dictionary()["atom_local_axes_matrix"][:, 0].reshape(3, 3)
    numpy.testing.assert_allclose(numpy.abs(matrix), numpy.eye(3), atol=1e-12)
    with pytest.raises(CrysPyError, match="has multipoles already"):
        crystal.add_multipole_form_factor("Cu")
    with pytest.raises(CrysPyError, match="no atom 'Cu1'"):
        crystal.add_multipole_form_factor("Cu1")


def test_mixed_atoms_form_factor_shape():
    """An atom with multipoles gets a form factor per operation, the other
    keeps its spherical one."""
    obj = cryspy.str_to_globaln(RCIF)
    crystal = obj.crystal_cu
    crystal.add_multipole_form_factor("Cu", coefficients="P20 0.3")
    dict_crystal = crystal.get_dictionary()
    index_hkl = numpy.array([[1., 0., 0.], [0., 0., 1.]]).transpose()
    rotation = numpy.stack([numpy.eye(3).flatten(), numpy.diag([0., 0., 0.]).flatten()], axis=1)
    rotation[:, 1] = numpy.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]).flatten()
    form_factor = numpy.full((2, 2), 0.7)
    result = calc_magnetic_form_factor_with_multipoles_by_dictionary(
        dict_crystal, index_hkl, rotation, numpy.array(["La", "Cu"]), form_factor)
    assert result.shape == (2, 2, 2)
    numpy.testing.assert_allclose(result[:, :, 0], 0.7)
    # P20 is anisotropic: the form factor along c differs from that along a*
    assert abs(result[0, 0, 1] - result[1, 0, 1]) > 1e-3


def test_assistant_finds_and_calls_the_method():
    from cryspy.procedure_assistant.methods import find_matching_methods
    obj = cryspy.str_to_globaln(RCIF)
    l_method = find_matching_methods(
        obj, "Add description of magnetic form factor through multipoles for Cu atom")
    assert l_method[0]["name"] == "add_multipole_form_factor"
    result = cryspy.call_assistant_tool(obj, "call_method", {
        "method": "add_multipole_form_factor", "arguments": {"atom_label": "Cu"}})
    assert "orbital 3d" in result["text"]
