"""Rows of loops and single items of any category added and removed by AI
assistants (add_item, remove_item): e.g. the multipoles of an atom."""
import pytest

import cryspy

RCIF = """global_

data_fe
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt "P 1"

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
  Fe  Fe3+   0.0  0.0  0.0  1.0  Biso  0.3
  D1  Dummy  0.0  0.0  0.5  1.0  Biso  0.0
  D2  Dummy  0.5  0.0  0.0  1.0  Biso  0.0
"""

MULTIPOLE_FE = {"atom_label": "Fe", "atom_type_symbol": "Fe3+", "orbital": "3d"}
# local axes of the atom, required by atom_rho_multipole (added with it)
AXES = {"atom_local_axes.atom0": "D1", "atom_local_axes.ax1": "Z",
        "atom_local_axes.atom2": "D2", "atom_local_axes.ax2": "X"}
AXES_D1 = {"atom_local_axes.atom0": "Fe", "atom_local_axes.ax1": "Z",
           "atom_local_axes.atom2": "D2", "atom_local_axes.ax2": "X"}


@pytest.fixture
def crystal():
    return cryspy.str_to_globaln(RCIF)


def call(obj, name, arguments):
    return cryspy.call_assistant_tool(obj, name, arguments)


def test_multipoles_of_an_atom_are_added_with_new_loops(crystal):
    result = call(crystal, "add_item", {"category": "atom_local_axes", "values": {
        "atom_label": "Fe", "atom0": "D1", "ax1": "Z", "atom1": "Fe", "atom2": "D2", "ax2": "X"}})
    assert not result["flag_error"] and result["flag_object_changed"], result["text"]
    # values written as JSON text, a coefficient by its CIF name
    result = call(crystal, "add_item", {"category": "atom_rho_multipole",
                                        "values": '{"atom_label": "Fe", "atom_type_symbol": "Fe3+", '
                                                  '"orbital": "3d", "coeff_P2-2": 0.1}'})
    assert "crystal_fe.atom_rho_multipole[Fe] is added:" in result["text"]
    # the magnetic ion needed by the multipoles is added with them, from atom_site
    assert "crystal_fe.atom_site_scat[Fe] is added: label Fe, type_symbol Fe3+" in result["text"]
    assert "atom_local_axes of Fe: there already, kept." in result["text"]
    fe = crystal.crystal_fe
    assert fe.atom_rho_multipole.items[0].coeff_P2m2 == pytest.approx(0.1)
    reloaded = cryspy.str_to_globaln(crystal.to_cif())
    assert [item.atom0 for item in reloaded.crystal_fe.atom_local_axes.items] == ["D1"]
    assert reloaded.crystal_fe.atom_rho_multipole.items[0].orbital == "3d"


def test_mandatory_values_are_asked_and_rows_keep_the_columns(crystal):
    result = call(crystal, "add_item", {"category": "atom_rho_multipole", "atom_label": "Fe"})
    # the ion comes from atom_site (through the new row of atom_site_scat),
    # the orbital from the knowledge base; the local axes are asked
    d_missing = {item["name"]: item for item in result["missing"]}
    assert list(d_missing) == ["atom_local_axes.atom0", "atom_local_axes.ax1",
                               "atom_local_axes.atom2", "atom_local_axes.ax2"]
    assert d_missing["atom_local_axes.atom0"]["proposal"] == "Fe: D1"
    assert "To add: atom_site_scat of Fe, atom_local_axes of Fe, atom_rho_multipole of Fe." in result["text"]
    call(crystal, "add_item", {"category": "atom_rho_multipole",
                               "values": dict(MULTIPOLE_FE, coeff_P20=0.2, **AXES)})
    call(crystal, "add_item", {"category": "atom_rho_multipole",
                               "values": dict(MULTIPOLE_FE, atom_label="D1", coeff_P22=0.3, **AXES_D1)})
    fe, d1 = crystal.crystal_fe.atom_rho_multipole.items
    # a coefficient of one row is 0 in the other one, not copied
    assert (fe.coeff_P22, d1.coeff_P20) == (0., 0.)
    reloaded = cryspy.str_to_globaln(crystal.to_cif()).crystal_fe.atom_rho_multipole.items
    assert [item.coeff_P22 for item in reloaded] == pytest.approx([0., 0.3])


@pytest.mark.parametrize("values, error", [
    (dict(MULTIPOLE_FE, atom_label="Co"), "no Co in atom_site.label"),
    (dict(MULTIPOLE_FE, colour=1), "atom_rho_multipole has no attribute colour"),
    (dict(MULTIPOLE_FE, kappa="big"), "'big' is not float"),
])
def test_wrong_values_are_explained(crystal, values, error):
    result = call(crystal, "add_item", {"category": "atom_rho_multipole", "values": values})
    assert result["flag_error"] and error in result["text"]


def test_existing_item_and_unknown_category_are_errors(crystal):
    call(crystal, "add_item", {"category": "atom_rho_multipole", "values": dict(MULTIPOLE_FE, **AXES)})
    result = call(crystal, "add_item", {"category": "atom_rho_multipole", "values": MULTIPOLE_FE})
    assert "Fe: atom_rho_multipole is there already, nothing is added." in result["text"]
    assert "crystal_fe.cell exists already" in call(crystal, "add_item", {"category": "cell"})["text"]
    assert "there is no category 'nothing'" in call(crystal, "add_item", {"category": "nothing"})["text"]
    result = call(crystal, "add_item", {"category": "extinction"})
    assert "no data block of the object can have extinction" in result["text"]


def test_data_block_is_asked_when_several_can_have_the_item(crystal):
    crystal.add_items([cryspy.str_to_globaln(RCIF.replace("data_fe", "data_co")).crystal_co])
    result = call(crystal, "add_item", {"category": "atom_rho_multipole", "values": MULTIPOLE_FE})
    assert result["missing"][0]["enum"] == ["crystal_fe", "crystal_co"]
    result = call(crystal, "add_item", {"category": "atom_rho_multipole", "block": "co",
                                        "values": dict(MULTIPOLE_FE, **AXES)})
    assert "crystal_co.atom_rho_multipole[Fe] is added" in result["text"]


def test_remove_rows_then_the_empty_loop(crystal):
    for label, axes in (("Fe", AXES), ("D1", AXES_D1)):
        call(crystal, "add_item", {"category": "atom_rho_multipole",
                                   "values": dict(MULTIPOLE_FE, atom_label=label, **axes)})
    result = call(crystal, "remove_item", {"path": "fe.atom_rho_multipole[D1]"})
    assert result["text"] == "crystal_fe.atom_rho_multipole[D1] is removed."
    assert len(crystal.crystal_fe.atom_rho_multipole.items) == 1
    result = call(crystal, "remove_item", {"path": "fe.atom_rho_multipole[Fe]"})
    assert result["text"] == "crystal_fe.atom_rho_multipole is removed."
    assert not crystal.crystal_fe.is_attribute("atom_rho_multipole")
    assert "cell is mandatory in crystal_fe" in call(crystal, "remove_item", {"path": "fe.cell"})["text"]


def test_help_on_multipoles_names_add_item():
    assert 'add_item {"category": "atom_rho_multipole"' in cryspy.get_help("multipoles")
