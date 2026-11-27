"""Conversion between single crystal FullProf files (Cry = 1, .int) and
rcif, checked against the chi2 of FullProf and by round trips."""
import os
import warnings

import pytest

import cryspy
from cryspy.H_functions_global.fullprof import (
    FullProfFormatError, load_fullprof, save_fullprof_files)

DIR_TESTS = os.path.dirname(os.path.abspath(__file__))
DIR_ROOT = os.path.dirname(os.path.dirname(DIR_TESTS))
FILE_HOTE3 = os.path.join(DIR_TESTS, "fullprof_single_crystal_hote3", "full.pcr")
FILE_CS3CR2BR9 = os.path.join(
    DIR_ROOT, "examples", "rhochi_unpolarized_single_crystal_k_Cs3Cr2Br9",
    "rhochi_unpolarized_single_crystal_k_Cs3Cr2Br9.rcif")
FILE_TBTE3 = os.path.join(DIR_TESTS, "twin_tbte3", "main.rcif")


def calc_chi_sq(rcif_object):
    result = cryspy.rhochi_no_refinement(rcif_object)
    return result["chi_sq"], result["n_point"]


def test_hote3_magnetic_as_fullprof():
    """HoTe3, k = (1/2, 1/2, 0), basis functions of Gamma1 (Jbt = 1, Isy = -2)
    and the extinction model 1: FullProf gives chi2 = 2170 for 116
    satellites. 0.5 %: CrysPy applies the extinction to each spin channel,
    which differs from FullProf by up to 1.4 % for l != 0."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rcif_object = load_fullprof(FILE_HOTE3)
    chi_sq, n_point = calc_chi_sq(rcif_object)
    assert n_point == 116
    assert chi_sq / (n_point - 1) == pytest.approx(2170., rel=5e-3)


@pytest.mark.parametrize("file_name", [FILE_HOTE3, FILE_CS3CR2BR9])
def test_magnetic_round_trip(file_name, tmp_path):
    """rcif -> FullProf -> rcif keeps chi_sq of the satellites."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if file_name.endswith(".rcif"):
            rcif_object = cryspy.load_file(file_name)
        else:
            rcif_object = load_fullprof(file_name)
        l_file = save_fullprof_files(rcif_object, str(tmp_path / "model.pcr"))
        assert os.path.isfile(str(tmp_path / "model.int"))
        rcif_object_2 = load_fullprof(l_file[0])
    chi_sq, n_point = calc_chi_sq(rcif_object)
    chi_sq_2, n_point_2 = calc_chi_sq(rcif_object_2)
    assert n_point_2 == n_point
    # the .pcr file keeps 4 digits of the moments
    assert chi_sq_2 == pytest.approx(chi_sq, rel=1e-3)


def test_nuclear_round_trip(tmp_path):
    """Nuclear structure (TbTe3 without twins and extinction): the same chi_sq
    and the refined scale."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rcif_object = cryspy.load_file(FILE_TBTE3)
        diffrn = [item for item in rcif_object.items if isinstance(item, cryspy.Diffrn)][0]
        diffrn.items[:] = [item for item in diffrn.items
                           if type(item).__name__ not in ("JanaTwinL", "Extinction")]
        l_file = save_fullprof_files(rcif_object, str(tmp_path / "model.pcr"))
        rcif_object_2 = load_fullprof(l_file[0])
    chi_sq, n_point = calc_chi_sq(rcif_object)
    chi_sq_2, n_point_2 = calc_chi_sq(rcif_object_2)
    assert n_point_2 == n_point
    assert chi_sq_2 == pytest.approx(chi_sq, rel=1e-3)
    assert len(rcif_object_2.get_variable_names()) == len(rcif_object.get_variable_names())


def test_flipping_ratios_are_refused(tmp_path):
    """Flipping ratios (itypdata = 2) are not converted yet."""
    text = open(FILE_HOTE3).read()
    with open(FILE_HOTE3.replace(".pcr", ".int")) as fid:
        l_line = fid.read().split("\n")
    l_line[2] = l_line[2].replace(" 0 ", " 2 ", 1)
    (tmp_path / "full.pcr").write_text(text)
    (tmp_path / "full.int").write_text("\n".join(l_line))
    with pytest.raises(FullProfFormatError):
        load_fullprof(str(tmp_path / "full.pcr"))
