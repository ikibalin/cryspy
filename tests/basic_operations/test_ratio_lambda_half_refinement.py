"""Refinement of the lambda/2 contamination (_setup_ratio_lambda/2)."""
import os

import numpy
import pytest

import cryspy

DIR_FUNCTIONAL = os.path.join(os.path.dirname(os.path.dirname(__file__)), "functional")
F_TBTE3 = os.path.join(DIR_FUNCTIONAL, "twin_tbte3", "main.rcif")


def test_refined_ratio_lambda_half_is_written_back():
    """Intensities simulated with ratio 0.03 give back 0.03 from a start of
    0.01: the refined value and its sigma reach the setup item."""
    obj = cryspy.load_file(F_TBTE3)
    diffrn = obj.diffrn_t1p8nuc
    diffrn.setup.ratio_lambdaover2 = 0.03
    cryspy.rhochi_no_refinement(obj)
    refln = diffrn.diffrn_refln
    intensity = numpy.array(refln.intensity_calc)
    refln.numpy_intensity = intensity
    refln.numpy_intensity_sigma = numpy.maximum(1e-3 * intensity, 1e-3)
    refln.numpy_to_items()
    obj.fix_variables()
    diffrn.setup.ratio_lambdaover2 = 0.01
    diffrn.setup.ratio_lambdaover2_refinement = True
    cryspy.rhochi_rietveld_refinement(obj)
    # rounding of intensity_calc in the items limits the agreement
    assert diffrn.setup.ratio_lambdaover2 == pytest.approx(0.03, abs=1e-4)
    assert diffrn.setup.ratio_lambdaover2_sigma > 0.
