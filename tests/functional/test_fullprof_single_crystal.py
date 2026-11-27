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


def write_hote3_extinction(tmp_path, line_value: str, line_code: str) -> str:
    """HoTe3 files with another extinction block of the .pcr file."""
    l_line = open(FILE_HOTE3).read().split("\n")
    index = [i for i, line in enumerate(l_line) if "Ext-Model" in line][0]
    l_line[index + 1], l_line[index + 2] = line_value, line_code
    (tmp_path / "full.pcr").write_text("\n".join(l_line))
    (tmp_path / "full.int").write_text(open(FILE_HOTE3.replace(".pcr", ".int")).read())
    return str(tmp_path / "full.pcr")


def calc_isotropic_q(ext_1: float):
    """q_hh, q_kk, q_ll of qani (and Ext1..Ext3 of the model 4) for the model 1
    of the HoTe3 cell: q_ii = Ext1/2 g*_ii."""
    from cryspy.H_functions_global.fullprof_single_crystal import calc_metric_reciprocal
    g_r = calc_metric_reciprocal([4.2645, 25.161797, 4.2799, 90., 90., 90.])
    return [ext_1 / 2. * g_r[i, i] for i in range(3)]


def test_extinction_model_4_as_model_1(tmp_path):
    """Ext-Model 4 (Ext1..Ext6 = q_hh ... q_kl of qani) with the values of the
    model 1 gives the same chi_sq, and it is written back as the model 1."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chi_sq_1, _ = calc_chi_sq(load_fullprof(FILE_HOTE3))
        file_pcr = write_hote3_extinction(
            tmp_path, "".join(f"{value:12.8f}" for value in calc_isotropic_q(0.5087))
            + "   0.0   0.0   0.0   0.0   4", "   0.00" * 7)
        rcif_object = load_fullprof(file_pcr)
        l_file = save_fullprof_files(rcif_object, str(tmp_path / "out" / "model.pcr"))
    # the cell of the .pcr file has more digits than this test
    assert calc_chi_sq(rcif_object)[0] == pytest.approx(chi_sq_1, rel=1e-4)
    assert open(l_file[0]).read().split("Ext-Model")[1].split("\n")[1].split()[-1] == "1"


def test_extinction_model_4(tmp_path):
    """Anisotropic Ext-Model 4: the refined coefficients with one code of
    FullProf are constrained, the model 4 is written back with the same
    chi_sq."""
    l_q = calc_isotropic_q(0.5087)
    l_q[2] *= 1.5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        file_pcr = write_hote3_extinction(
            tmp_path, "".join(f"{value:12.8f}" for value in l_q) + "   0.0   0.0   0.0   0.0   4",
            "  31.00  31.00  41.00   0.00   0.00   0.00   0.00")
        rcif_object = load_fullprof(file_pcr)
        chi_sq_4, _ = calc_chi_sq(rcif_object)
        extinction = [item for item in rcif_object.items
                      if isinstance(item, cryspy.Diffrn)][0].extinction
        assert extinction.q_ll == pytest.approx(l_q[2])
        assert extinction.q_hh_refinement and extinction.q_ll_refinement
        l_expression = [constraint.expression for item in rcif_object.items
                        if type(item).__name__ == "ConstraintL" for constraint in item.items]
        assert len(l_expression) == 1 and l_expression[0].startswith("ext_kk = 1*ext_hh")
        l_file = save_fullprof_files(rcif_object, str(tmp_path / "out" / "model.pcr"))
        text = open(l_file[0]).read()
        assert text.split("Ext-Model")[1].split("\n")[1].split()[-1] == "4"
        chi_sq_back, _ = calc_chi_sq(load_fullprof(l_file[0]))
    assert chi_sq_back == pytest.approx(chi_sq_4, rel=1e-4)
