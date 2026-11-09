"""Structures from COD (cryspy.H_functions_global.cod), without network."""
import json

import numpy
import pytest

import cryspy
from cryspy.H_functions_global import cod

CIF_COD = """#------------------------------------------------------------------------------
data_9008460
_chemical_formula_sum            Al
_symmetry_space_group_name_H-M   'F m -3 m'
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
_cell_length_a                   4.04958
_cell_length_b                   4.04958
_cell_length_c                   4.04958
_cod_database_code               9008460
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Al1 0.00000 0.00000 0.00000
"""


@pytest.mark.parametrize("formula, hill", [
    ("Al", "Al"), ("TiO2", "O2 Ti"), ("C2H6O", "C2 H6 O"), ("Te3Tb", "Tb Te3")])
def test_hill_formula(formula, hill):
    assert cod.form_hill_formula(formula) == hill


def test_cif_to_crystal():
    crystal = cod.cif_to_crystal(CIF_COD, data_name="al")
    atom = crystal.atom_site.items[0]
    assert atom.type_symbol == "Al"
    assert atom.adp_type == "Biso"
    # F(111) of fcc Al: 4 b(Al), b(Al) = 0.3449 10^-12 cm
    f_111 = crystal.calc_f_nucl(numpy.array([[1], [1], [1]]))
    assert abs(f_111[0]) == pytest.approx(4 * 0.3449, abs=1e-3)


def test_add_crystal_from_cod(monkeypatch):
    l_url = []

    def read_url(url, data=None):
        l_url.append(url)
        if "result?" in url:
            return json.dumps([{"file": "9008460"}])
        return CIF_COD

    monkeypatch.setattr(cod, "read_url", read_url)
    obj = cryspy.str_to_globaln("global_")
    result = cryspy.call_assistant_tool(
        obj, "add_crystal_from_cod", {"formula": "Al"})
    assert not result["flag_error"], result["text"]
    assert result["flag_object_changed"]
    assert l_url[-1].endswith("/9008460.cif")
    cryspy.call_assistant_tool(obj, "add_crystal_from_cod", {"cod_id": "9008460"})
    l_name = [item.data_name for item in obj.items]
    assert l_name == ["al", "cod_9008460"]
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    assert len(obj_2.items) == 2


def test_cod_errors_are_text(monkeypatch):
    def read_url(url, data=None):
        raise OSError("no network")

    monkeypatch.setattr(cod, "read_url", read_url)
    obj = cryspy.str_to_globaln("global_")
    result = cryspy.call_assistant_tool(
        obj, "add_crystal_from_cod", {"formula": "Al"})
    assert result["flag_error"]
    assert "COD cannot be reached" in result["text"]


def test_cif_without_space_group_and_with_u_iso():
    cif_text = """data_x
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
O1 O2- 0.1 0.2 0.3 0.01
"""
    crystal = cod.cif_to_crystal(cif_text)
    atom = crystal.atom_site.items[0]
    assert atom.adp_type == "Uiso"
    assert atom.type_symbol == "O2-"
    assert len(crystal.space_group.full_space_group_symop.items) == 1
