"""Objects of CIF categories unknown to CrysPy (generic ItemN and LoopN,
e.g. _publ_author_name of a COD file) are copied and pickled.

copy.deepcopy and pickle make the object without __init__ and ask it for
__setstate__; __getattr__ looked for ATTR_NAMES / ITEM_CLASS, which such a
generic object keeps in its own __dict__, so it called itself until
RecursionError. Crystal.report_magnetic_models copies the crystal: in a
thread of CrysPy Editor on macOS (512 KB stack) the editor closed with
'bus error' for La2CuO4 from COD."""
import copy
import pickle

import cryspy
from cryspy.B_parent_classes.cl_1_item import is_generic

CIF = """data_la2cuo4
loop_
_publ_author_name
'Tuilier, M H'
'Chevalier, B'
_journal_year 1992
_space_group_name_H-M_alt 'C m c e'
_cell_length_a 5.352
_cell_length_b 5.4
_cell_length_c 13.157
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_occupancy
La1 La 0. 0. 0.3617 Biso 0.5 1.0
Cu1 Cu 0. 0. 0. Biso 0.3 1.0
O1 O 0.25 0.25 0.0 Biso 0.5 1.0
O2 O 0. 0. 0.182 Biso 0.5 1.0
"""


def take_crystal():
    obj = cryspy.str_to_globaln(CIF)
    crystal = obj.crystal_la2cuo4
    generic = [item for item in crystal.items if is_generic(item)]
    assert generic, "the test needs an unknown CIF category"
    return obj, crystal


def test_deepcopy_with_unknown_categories():
    obj, crystal = take_crystal()
    for item in (obj, crystal) + tuple(crystal.items):
        assert copy.deepcopy(item).to_cif() == item.to_cif()


def test_pickle_with_unknown_categories():
    obj, crystal = take_crystal()
    assert pickle.loads(pickle.dumps(crystal)).to_cif() == crystal.to_cif()


def test_report_magnetic_models_of_a_cod_file():
    _, crystal = take_crystal()
    report = crystal.report_magnetic_models(0., 0., 0., "Cu1")
    assert "| 1 | Γ1 |" in report and "not possible" not in report
