"""The unit cell of a paramagnetic crystal is refined against flipping ratios
(the derivative by eq_ccs failed with KeyError 'eq_css')."""
import os

import pytest

import cryspy

F_EXAMPLE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "examples",
    "rhochi_single_susceptibility_Ho2Ti2O7", "rhochi_single_susceptibility_Ho2Ti2O7.rcif")


def test_refine_unit_cell_by_flip_ratio():
    cryspy_object = cryspy.load_file(F_EXAMPLE)
    crystal = [item for item in cryspy_object.items if isinstance(item, cryspy.Crystal)][0]
    crystal.cell.length_a_refinement = True
    chi_sq_start = cryspy.rhochi_no_refinement(cryspy_object)["chi_sq"]
    result = cryspy.rhochi_rietveld_refinement(cryspy_object)
    assert result["chi_sq"] <= chi_sq_start * (1. + 1e-6)
    assert crystal.cell.length_a == pytest.approx(crystal.cell.length_b)
