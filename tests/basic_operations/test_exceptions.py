"""CrysPy errors: own classes with the former built-in classes as bases."""
import os
import warnings

import pytest

import cryspy

CELL = """
_cell_length_a {a:}
_cell_length_b 5.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
"""

ATOM_SITE = """
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
Fe Fe 0.0 {y:} 0.25 {adp:}
"""


def test_hierarchy():
    assert issubclass(cryspy.RcifFormatError, cryspy.CrysPyError)
    assert issubclass(cryspy.RcifFormatError, ValueError)
    assert issubclass(cryspy.RcifFileNotFoundError, cryspy.CrysPyError)
    assert issubclass(cryspy.RcifFileNotFoundError, FileNotFoundError)
    assert issubclass(cryspy.RefinementError, cryspy.CrysPyError)
    assert issubclass(cryspy.RefinementError, AttributeError)
    assert issubclass(cryspy.CrysPyWarning, UserWarning)


def test_file_not_found(tmp_path):
    f_name = os.path.join(str(tmp_path), "absent.rcif")
    with pytest.raises(cryspy.RcifFileNotFoundError, match="absent.rcif"):
        cryspy.file_to_globaln(f_name)
    with pytest.raises(FileNotFoundError):
        cryspy.Cell.from_cif_file(f_name)


def test_broken_text():
    with pytest.raises(cryspy.RcifFormatError):
        cryspy.str_to_globaln("this is not a cif")


def test_not_a_number_in_item():
    with pytest.raises(cryspy.RcifFormatError, match="_cell_length_a"):
        cryspy.str_to_globaln("data_x" + CELL.format(a="abc"))


def test_not_a_number_in_loop():
    with pytest.raises(cryspy.RcifFormatError, match="_atom_site_fract_y"):
        cryspy.str_to_globaln("data_x" + ATOM_SITE.format(y="x", adp="Uiso"))


def test_not_an_integer():
    with pytest.raises(cryspy.RcifFormatError, match="index_h"):
        cryspy.str_to_globaln(
            "data_x\nloop_\n_diffrn_refln_index_h\n_diffrn_refln_index_k\n"
            "_diffrn_refln_index_l\n_diffrn_refln_fr\n_diffrn_refln_fr_sigma\n"
            "0.5 0 0 1.1 0.1\n")


def test_unknown_value_and_question_mark_are_accepted():
    with warnings.catch_warnings():
        warnings.simplefilter("error", cryspy.CrysPyWarning)
        for value in (".", "?"):
            rcif = cryspy.str_to_globaln("data_x" + ATOM_SITE.format(y=value, adp="Uiso"))
            assert rcif.data_x.atom_site["Fe"].fract_z == 0.25


def test_value_outside_allowed_list_warns():
    with pytest.warns(cryspy.CrysPyWarning, match="adp_type"):
        rcif = cryspy.str_to_globaln("data_x" + ATOM_SITE.format(y="0.5", adp="Bxyz"))
    assert not rcif.data_x.atom_site["Fe"].is_attribute("adp_type")


def test_refinement_without_crystal():
    rcif = cryspy.str_to_globaln("data_x" + CELL.format(a="4.0"))
    with pytest.raises(cryspy.RefinementError, match="Crystal"):
        cryspy.rhochi_rietveld_refinement(rcif)
    with pytest.raises(AttributeError):
        cryspy.rhochi_rietveld_refinement(rcif)
