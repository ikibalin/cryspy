"""Matching mode of 1D powder patterns (pd_matching)."""
import os

import numpy
import pytest

import cryspy
from cryspy.C_item_loop_classes.cl_1_pd_peak import PdPeakL

DIR = os.path.dirname(__file__)
F_PBSO4 = os.path.join(DIR, "PbSO4_unpol_powder_test", "main.rcif")
F_HO2TI2O7 = os.path.join(DIR, "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif")


def take_block(rcif, class_name: str):
    return [block for block in rcif.items if type(block).__name__ == class_name][0]


def take_peaks(experiment, label: str):
    return [item for item in experiment.items
            if isinstance(item, PdPeakL) and item.loop_name == label][0]


def load_pbso4(ttheta_max: float = None):
    rcif = cryspy.file_to_globaln(F_PBSO4)
    experiment = take_block(rcif, "Pd")
    if ttheta_max is not None:
        experiment.range.ttheta_max = ttheta_max
    return rcif, experiment, take_block(rcif, "Crystal")


def test_fixed_intensities_reproduce_rietveld():
    rcif, experiment, _ = load_pbso4()
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    experiment.switch_matching_mode("fixed_intensities")
    chi_sq_fixed = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    assert chi_sq_fixed == pytest.approx(chi_sq_rietveld, rel=1e-6)


def test_refine_intensities_unpolarized():
    rcif, experiment, crystal = load_pbso4()
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    experiment.switch_matching_mode("refine_intensities")
    # the atoms are not used: the result does not depend on them
    crystal.atom_site.items[0].fract_x = 0.3
    chi_sq_matching = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    assert chi_sq_matching < 0.6 * chi_sq_rietveld
    peaks = take_peaks(experiment, "pbso4")
    intensity_plus = peaks.numpy_intensity_plus
    assert numpy.all(intensity_plus >= 0.)
    assert numpy.any(intensity_plus > 0.)
    # unpolarized neutrons: the halves up and down are equal
    assert numpy.allclose(intensity_plus, peaks.numpy_intensity_minus)

    # the fitted intensities kept fixed give the same pattern
    rcif_2 = cryspy.str_to_globaln(rcif.to_cif())
    experiment_2 = take_block(rcif_2, "Pd")
    assert experiment_2.pd_matching.mode == "refine_intensities"
    experiment_2.switch_matching_mode("fixed_intensities")
    chi_sq_fixed = cryspy.rhochi_no_refinement(rcif_2)["chi_sq"]
    # pd_peak keeps two decimals of the intensities
    assert chi_sq_fixed == pytest.approx(chi_sq_matching, rel=1e-4)

    experiment.switch_matching_mode("off")
    assert cryspy.rhochi_no_refinement(rcif)["chi_sq"] > chi_sq_rietveld


def test_refine_intensities_and_zero_shift():
    rcif, experiment, _ = load_pbso4(ttheta_max=60.)
    rcif.fix_variables()
    offset_true = experiment.setup.offset_ttheta
    experiment.switch_matching_mode("refine_intensities")
    experiment.setup.offset_ttheta = offset_true + 0.05
    chi_sq_start = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    experiment.setup.offset_ttheta_refinement = True
    result = cryspy.rhochi_rietveld_refinement(rcif)
    # the zero shift of the Rietveld refinement is found again with free
    # intensities
    assert experiment.setup.offset_ttheta == pytest.approx(offset_true, abs=0.02)
    assert result["chi_sq"] < 0.3 * chi_sq_start


def refine_length_a(mode: str):
    rcif, experiment, crystal = load_pbso4(ttheta_max=60.)
    rcif.fix_variables()
    cryspy.rhochi_no_refinement(rcif)
    experiment.switch_matching_mode(mode)
    intensity_plus = take_peaks(experiment, "pbso4").numpy_intensity_plus
    crystal.cell.length_a = crystal.cell.length_a * 1.001
    crystal.cell.length_a_refinement = True
    cryspy.rhochi_rietveld_refinement(rcif)
    return crystal.cell.length_a, intensity_plus, take_peaks(experiment, "pbso4").numpy_intensity_plus


def test_fixed_intensities_and_cell():
    length_a_rietveld = refine_length_a("off")[0]
    length_a, intensity_start, intensity_end = refine_length_a("fixed_intensities")
    # the intensities of the Rietveld calculation kept fixed: the same cell
    # (the structure factors depend on the cell only weakly)
    assert length_a == pytest.approx(length_a_rietveld, abs=2e-4)
    assert length_a != pytest.approx(8.479343 * 1.001, abs=1e-3)
    # the intensities are not changed
    assert numpy.allclose(intensity_end, intensity_start, atol=0.01)


def test_fixed_intensities_need_pd_peak():
    rcif, experiment, _ = load_pbso4()
    experiment.items[:] = [item for item in experiment.items if not isinstance(item, PdPeakL)]
    with pytest.raises(ValueError, match="pd_peak"):
        experiment.switch_matching_mode("fixed_intensities")
    with pytest.raises(ValueError, match="matching mode"):
        experiment.switch_matching_mode("le_bail")


def test_polarized_powder():
    rcif = cryspy.file_to_globaln(F_HO2TI2O7)
    experiment = take_block(rcif, "Pd")
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    experiment.switch_matching_mode("fixed_intensities")
    assert cryspy.rhochi_no_refinement(rcif)["chi_sq"] == pytest.approx(chi_sq_rietveld, rel=1e-6)

    experiment.switch_matching_mode("refine_intensities")
    result = cryspy.rhochi_no_refinement(rcif)
    assert result["chi_sq"] < 0.3 * chi_sq_rietveld
    peaks = take_peaks(experiment, "phase1")
    intensity_plus = peaks.numpy_intensity_plus
    intensity_minus = peaks.numpy_intensity_minus
    assert numpy.all(intensity_plus >= 0.) and numpy.all(intensity_minus >= 0.)
    # the field polarizes Ho: up and down differ
    assert not numpy.allclose(intensity_plus, intensity_minus)
    assert "# Matching mode" in experiment.report()


def test_recipe_profile_matching():
    """The assistant recipe: matching mode switched on by call_method, the
    zero shift, the cell and the widths refined with free intensities."""
    rcif, experiment, _ = load_pbso4(ttheta_max=40.)
    result = cryspy.run_recipe(rcif, "profile_matching")
    d_status = {name: status for name, status, _ in result["steps"]}
    assert d_status["start"] == "done", result["text"]
    assert experiment.pd_matching.mode == "refine_intensities"
    assert d_status["zero"] == "done" and d_status["cell"] == "done", result["text"]
    assert numpy.all(take_peaks(experiment, "pbso4").numpy_intensity_plus >= 0.)


def test_diagnose_fit_in_matching_mode():
    rcif, experiment, _ = load_pbso4(ttheta_max=40.)
    experiment.switch_matching_mode("refine_intensities")
    cryspy.call_assistant_tool(rcif, "calculate")
    text = cryspy.call_assistant_tool(rcif, "diagnose_fit")
    text = text["text"] if isinstance(text, dict) else text
    assert "MATCHING" in text
    assert "SCALE" not in text and "INTENSITY_SCATTER" not in text
