"""Conversion between FullProf files and rcif, checked against the
calculations of FullProf (.prf files) and by round trips."""
import os
import warnings

import numpy
import pytest

import cryspy
from cryspy.H_functions_global.fullprof import (
    FullProfFormatError, load_fullprof, read_prf_columns, save_fullprof_files)

DIR_TESTS = os.path.dirname(os.path.abspath(__file__))
DIR_ROOT = os.path.dirname(os.path.dirname(DIR_TESTS))
DIR_PBSO4 = os.path.join(DIR_TESTS, "PbSO4_unpol_powder_test")
DIR_TWO_PHASES = os.path.join(DIR_TESTS, "two_phases_simul_pd", "fp_nomag_2phases")
DIR_XRAYS = os.path.join(
    DIR_ROOT, "examples", "rhochi_unpolarized_powder_1d_xrays_PbSO4", "fullprof_xray")


def take_pd(rcif_object):
    return [item for item in rcif_object.items if isinstance(item, cryspy.Pd)][0]


def calc_r_net(rcif_object, file_prf, ttheta_min, ttheta_max) -> float:
    """R factor between the net calculated patterns (without background) of
    CrysPy and FullProf."""
    cryspy.rhochi_no_refinement(rcif_object)
    pd_proc = take_pd(rcif_object).pd_proc
    ttheta = numpy.array(pd_proc.ttheta, dtype=float)
    net = (numpy.array(pd_proc.intensity_plus_net, dtype=float)
           + numpy.array(pd_proc.intensity_minus_net, dtype=float))
    ttheta_fp, y_calc_fp, background_fp = read_prf_columns(file_prf)
    flag = (ttheta_fp >= max(ttheta_min, ttheta.min())) & (ttheta_fp <= min(ttheta_max, ttheta.max()))
    net_fp = y_calc_fp[flag] - background_fp[flag]
    net_cryspy = numpy.interp(ttheta_fp[flag], ttheta, net)
    return numpy.abs(net_fp - net_cryspy).sum() / numpy.abs(net_fp).sum()


def test_xray_simulation_as_fullprof():
    """X-ray simulation of FullProf (Job = 2): the converted model gives the
    same pattern. FullProf has no asymmetry here, so the comparison checks
    the structure factors, occupancies, Lorentz-polarization factor and
    the X-ray scale (2 x FullProf)."""
    file_pcr = os.path.join(DIR_XRAYS, "structure.pcr")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rcif_object = load_fullprof(file_pcr)
    file_prf = os.path.join(DIR_XRAYS, "structure.prf")
    for ttheta_min, ttheta_max in ((0., 40.), (40., 80.), (80., 180.)):
        # 5 %: the patterns differ by a scale of 2-3 % (X-ray form factors of the
        # ions), the rest (1-2 %) is the interpolation and the cut-off of the tails
        assert calc_r_net(rcif_object, file_prf, ttheta_min, ttheta_max) < 0.05


def test_refined_pbso4_neutrons():
    """The final .pcr file of FullProf (pbso4.new, chi2 4.0-4.2) gives the same
    fit in CrysPy after the refinement of the zero shift and the empirical
    asymmetry, whose convention differs (issue #50)."""
    file_pcr = os.path.join(DIR_PBSO4, "fp", "pbso4.new")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rcif_object = load_fullprof(file_pcr)
    rcif_object.fix_variables()
    pd = take_pd(rcif_object)
    pd.setup.offset_ttheta_refinement = True
    asymmetry = pd.pd_instr_reflex_asymmetry
    asymmetry.p1_refinement = True
    asymmetry.p2_refinement = True
    asymmetry.p3_refinement = True
    cryspy.rhochi_rietveld_refinement(rcif_object)
    result = cryspy.rhochi_no_refinement(rcif_object)
    chi_sq_per_point = result["chi_sq"] / result["n_point"]
    assert chi_sq_per_point < 4.5
    # FullProf: Zero = -0.088 degrees with the same sign convention
    assert pd.setup.offset_ttheta == pytest.approx(-0.088, abs=0.01)


def test_two_phases_with_betas():
    """Two phases, anisotropic displacements (betas), background of the
    .prf file (Nba = -1) and AsyLim."""
    file_pcr = os.path.join(DIR_TWO_PHASES, "full.pcr")
    with warnings.catch_warnings(record=True) as l_warning:
        warnings.simplefilter("always")
        rcif_object = load_fullprof(file_pcr)
    l_message = [str(warning.message) for warning in l_warning]
    assert any("AsyLim" in message for message in l_message)
    l_crystal = [item for item in rcif_object.items if isinstance(item, cryspy.Crystal)]
    assert len(l_crystal) == 2
    assert any(crystal.is_attribute("atom_site_aniso") for crystal in l_crystal)
    pd = take_pd(rcif_object)
    assert len(pd.phase.items) == 2
    assert len(pd.pd_background.items) > 0


@pytest.mark.parametrize("file_name", [
    os.path.join(DIR_PBSO4, "main.rcif"),
    os.path.join(DIR_PBSO4, "fp", "pbso4.pcr"),
    os.path.join(DIR_TWO_PHASES, "full.pcr")])
def test_round_trip(file_name, tmp_path):
    """rcif -> FullProf -> rcif keeps chi_sq and the refined parameters."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if file_name.endswith(".rcif"):
            rcif_object = cryspy.load_file(file_name)
        else:
            rcif_object = load_fullprof(file_name)
        l_file = save_fullprof_files(rcif_object, str(tmp_path / "model.pcr"))
        assert os.path.isfile(str(tmp_path / "model.dat"))
        rcif_object_2 = load_fullprof(l_file[0])
    result = cryspy.rhochi_no_refinement(rcif_object)
    result_2 = cryspy.rhochi_no_refinement(rcif_object_2)
    assert result_2["n_point"] == result["n_point"]
    # the .pcr file keeps 6 significant digits
    assert result_2["chi_sq"] == pytest.approx(result["chi_sq"], rel=1e-3)
    assert len(rcif_object_2.get_variable_names()) == len(rcif_object.get_variable_names())


def test_magnetic_phase_is_refused():
    """FullProf magnetic phases (Jbt = 10) are not converted."""
    file_pcr = os.path.join(
        DIR_TESTS, "IG_Fe_core_Mn_shell_powder1d_fp_test")
    l_file = [os.path.join(directory, name) for directory, _, l_name in os.walk(file_pcr)
              for name in l_name if name.endswith(".pcr")]
    assert len(l_file) > 0
    with pytest.raises(FullProfFormatError):
        load_fullprof(l_file[0])


def test_procedures_and_tools(tmp_path):
    """add_fullprof_files and save_fullprof_files as procedures of the
    object and as tools of AI assistants."""
    from cryspy.procedure_assistant.tools import call_assistant_tool
    rcif_object = cryspy.GlobalN.make_container((), (), "global_")
    text = call_assistant_tool(rcif_object, "add_fullprof_files", {
        "file_pcr": os.path.join(DIR_PBSO4, "fp", "pbso4.pcr")})
    assert "added" in text["text"], text
    assert any(isinstance(item, cryspy.Pd) for item in rcif_object.items)
    text = call_assistant_tool(rcif_object, "save_fullprof_files", {
        "file_pcr": str(tmp_path / "copy.pcr")})
    assert os.path.isfile(str(tmp_path / "copy.pcr")), text
    with pytest.raises(ValueError):
        cryspy.add_fullprof_files(rcif_object, os.path.join(DIR_PBSO4, "fp", "pbso4.pcr"))
