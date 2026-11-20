"""Reading a file gives the objects of pycifstar to the classes of CrysPy
(from_cif_data) instead of writing them to text and reading it again for
every class: faster, and the text values are kept as they are."""
import pytest

import cryspy
from cryspy.A_functions_base.function_2_space_group import \
    get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell

COD = """data_x
_journal_name_full 'Physica C'
_publ_section_abstract ?
_publ_section_title
; ' CeCuAl3 RT #13744 ThCr2Si2 type'
;
_publ_section_title_footnote
;
;
_cell_length_a 4.2
_cell_length_b 4.2
_cell_length_c 10.6
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'I 4/m m m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ce1 Ce 0 0 0
"""


def test_text_with_hash_is_kept():
    """'#' inside a text field is not a comment: before, the title was
    cut to 'CeCuAl3 RT' when the data were written and read again."""
    obj = cryspy.str_to_globaln(COD.split("_cell_length_a")[0])
    publ = obj.items[0].publ_section
    assert publ.title.strip() == "CeCuAl3 RT #13744 ThCr2Si2 type"
    # an empty text field is written as '.', not as a tag without a value
    text = obj.to_cif()
    assert "_publ_section_title_footnote .\n" in text


def test_from_cif_data_is_from_cif():
    from pycifstar import Data
    data = Data()
    data.take_from_string(COD)
    assert Cell.from_cif_data(data).to_cif() == Cell.from_cif(COD).to_cif()
    assert ItemN.from_cif_data(Data(items=data.items.items_with_prefix("_journal"))).name_full == \
        "Physica C"
    loop = LoopN.from_cif_data(Data(loops=list(data.loops)))
    assert loop.label == ["Ce1"]


def test_cached_space_group_tables_are_copies():
    tables = get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2(139, "1")
    tables[0].append("x,y,z")
    tables[2][0] = -1
    again = get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2(139, "1")
    assert again[0] != tables[0] and again[2][0] != -1


def test_writing_attribute_drops_internal_parameters():
    cell = Cell(length_a=5., length_b=5., length_c=5., angle_alpha=90.,
                angle_beta=90., angle_gamma=90.)
    cell.m_g
    cell.__dict__["file_input"] = "x.rcif"  # not a CIF attribute
    cell.length_a = 4.
    assert "m_g" not in cell.__dict__ and "file_input" not in cell.__dict__
    assert cell.m_g[0, 0] == pytest.approx(16.)
