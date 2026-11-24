"""Profile matching of 1D powder patterns (rhochi_profile_matching)."""
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
    rcif.fix_variables()
    return rcif, experiment, take_block(rcif, "Crystal")


def test_fixed_intensities_reproduce_rietveld():
    rcif, _, _ = load_pbso4()
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    chi_sq_fixed = cryspy.rhochi_profile_matching(rcif, fixed_intensities=True)["chi_sq"]
    assert chi_sq_fixed == pytest.approx(chi_sq_rietveld, rel=1e-6)


def test_free_intensities_unpolarized():
    rcif, experiment, crystal = load_pbso4()
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    # the atoms are not used: the result does not depend on them
    crystal.atom_site.items[0].fract_x = 0.3
    chi_sq_matching = cryspy.rhochi_profile_matching(rcif)["chi_sq"]
    assert chi_sq_matching < 0.6 * chi_sq_rietveld
    peaks = take_peaks(experiment, "pbso4")
    intensity_plus = peaks.numpy_intensity_plus
    assert numpy.all(intensity_plus >= 0.)
    assert numpy.any(intensity_plus > 0.)
    # unpolarized neutrons: the halves up and down are equal
    assert numpy.allclose(intensity_plus, peaks.numpy_intensity_minus)
    # the procedure writes nothing about itself to the experiment
    s_cif = rcif.to_cif()
    assert "matching" not in s_cif

    # the fitted intensities kept fixed give the same pattern after the rcif
    rcif_2 = cryspy.str_to_globaln(s_cif)
    chi_sq_fixed = cryspy.rhochi_profile_matching(rcif_2, fixed_intensities=True)["chi_sq"]
    # pd_peak keeps two decimals of the intensities
    assert chi_sq_fixed == pytest.approx(chi_sq_matching, rel=1e-4)

    # the next calculation is a Rietveld one again
    assert cryspy.rhochi_no_refinement(rcif)["chi_sq"] > chi_sq_rietveld


def test_free_intensities_and_zero_shift():
    rcif, experiment, _ = load_pbso4(ttheta_max=60.)
    offset_true = experiment.setup.offset_ttheta
    experiment.setup.offset_ttheta = offset_true + 0.05
    chi_sq_start = cryspy.rhochi_profile_matching(rcif)["chi_sq"]
    experiment.setup.offset_ttheta_refinement = True
    result = cryspy.rhochi_profile_matching(rcif)
    # the zero shift of the Rietveld refinement is found again with free
    # intensities, with its error
    assert experiment.setup.offset_ttheta == pytest.approx(offset_true, abs=0.02)
    assert experiment.setup.offset_ttheta_sigma > 0.
    assert result["chi_sq"] < 0.3 * chi_sq_start


def refine_length_a(fixed_intensities):
    """Cell length a refined from a shifted value by the Rietveld refinement
    (None) or the profile matching with the intensities of the Rietveld
    calculation kept fixed."""
    rcif, experiment, crystal = load_pbso4(ttheta_max=60.)
    cryspy.rhochi_no_refinement(rcif)
    intensity_plus = take_peaks(experiment, "pbso4").numpy_intensity_plus
    crystal.cell.length_a = crystal.cell.length_a * 1.001
    crystal.cell.length_a_refinement = True
    if fixed_intensities is None:
        cryspy.rhochi_rietveld_refinement(rcif)
    else:
        cryspy.rhochi_profile_matching(rcif, fixed_intensities=fixed_intensities)
    return crystal.cell.length_a, intensity_plus, take_peaks(experiment, "pbso4").numpy_intensity_plus


def test_fixed_intensities_and_cell():
    length_a_rietveld = refine_length_a(None)[0]
    length_a, intensity_start, intensity_end = refine_length_a(True)
    # the same cell (the structure factors depend on the cell only weakly)
    assert length_a == pytest.approx(length_a_rietveld, abs=2e-4)
    assert length_a != pytest.approx(8.479343 * 1.001, abs=1e-3)
    # the intensities are not changed
    assert numpy.allclose(intensity_end, intensity_start, atol=0.01)


def test_fixed_intensities_need_pd_peak():
    rcif, experiment, _ = load_pbso4()
    experiment.items[:] = [item for item in experiment.items if not isinstance(item, PdPeakL)]
    with pytest.raises(ValueError, match="pd_peak"):
        cryspy.rhochi_profile_matching(rcif, fixed_intensities=True)


def test_polarized_powder():
    rcif = cryspy.file_to_globaln(F_HO2TI2O7)
    rcif.fix_variables()
    experiment = take_block(rcif, "Pd")
    chi_sq_rietveld = cryspy.rhochi_no_refinement(rcif)["chi_sq"]
    chi_sq_fixed = cryspy.rhochi_profile_matching(rcif, fixed_intensities=True)["chi_sq"]
    assert chi_sq_fixed == pytest.approx(chi_sq_rietveld, rel=1e-6)

    result = cryspy.rhochi_profile_matching(rcif)
    assert result["chi_sq"] < 0.3 * chi_sq_rietveld
    peaks = take_peaks(experiment, "phase1")
    intensity_plus = peaks.numpy_intensity_plus
    intensity_minus = peaks.numpy_intensity_minus
    assert numpy.all(intensity_plus >= 0.) and numpy.all(intensity_minus >= 0.)
    # the field polarizes Ho: up and down differ
    assert not numpy.allclose(intensity_plus, intensity_minus)


def test_procedure_for_menus_and_assistants():
    description = cryspy.describe_procedure(cryspy.rhochi_profile_matching)
    l_name = [d_param["name"] for d_param in description["parameters"]]
    assert "fixed_intensities" in l_name and "callback" not in l_name


def test_assistant_tool_and_recipe():
    """The tool profile_matching (undo works as after refine) and the recipe:
    the zero shift and the cell refined with free intensities."""
    rcif, experiment, _ = load_pbso4(ttheta_max=40.)
    result = cryspy.call_assistant_tool(rcif, "profile_matching")
    assert not result["flag_error"], result["text"]
    assert "pd_peak" in result["text"]

    result = cryspy.run_recipe(rcif, "profile_matching")
    d_status = {name: status for name, status, _ in result["steps"]}
    assert d_status["start"] == "done", result["text"]
    assert d_status["zero"] in ("done", "undone") and d_status["cell"] in ("done", "undone"), result["text"]
    assert numpy.all(take_peaks(experiment, "pbso4").numpy_intensity_plus >= 0.)

    rcif, _, _ = load_pbso4(ttheta_max=40.)
    result = cryspy.run_recipe(rcif, "refine_with_fixed_intensities")
    d_status = {name: status for name, status, _ in result["steps"]}
    assert d_status["start"] == "done", result["text"]
