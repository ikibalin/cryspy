"""Sigma of a symmetry-constrained susceptibility after a refinement.

On the cubic site of Ho2Ti2O7 the constraint chi_11 = chi_22 = chi_33 is a
projector with weights 1/3: the refined chi_11 enters each component of the
tensor with the weight 1/3. Its sigma was written to chi_11 unscaled (three
times too large) and nothing to chi_22, chi_33.
"""
import os

import numpy
import pytest

import cryspy

F_NAME = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif")
COMPONENTS = ("11", "22", "33", "12", "13", "23")


def refine(component_1: str, component_2: str):
    obj = cryspy.file_to_globaln(F_NAME)
    obj.fix_variables()
    crystal = [item for item in obj.items if isinstance(item, cryspy.Crystal)][0]
    chi = crystal.atom_site_susceptibility.items[0]
    for name in (component_1, component_2):
        setattr(chi, f"chi_{name}_refinement", True)
    result = cryspy.rhochi_rietveld_refinement(obj)
    # the object reproduces the chi_sq of the minimizer
    assert cryspy.rhochi_no_refinement(obj)["chi_sq"] == pytest.approx(
        result["chi_sq"], rel=1e-9)
    return (numpy.array([float(getattr(chi, f"chi_{c}")) for c in COMPONENTS]),
            numpy.array([float(getattr(chi, f"chi_{c}_sigma")) for c in COMPONENTS]))


def test_sigma_goes_to_all_equivalent_components():
    value, sigma = refine("11", "12")
    assert numpy.all(sigma > 0.)
    numpy.testing.assert_allclose(sigma[:3], sigma[0])
    numpy.testing.assert_allclose(sigma[3:], sigma[3])
    numpy.testing.assert_allclose(value[:3], value[0])

    # refining another component of the same group gives the same answer
    value_2, sigma_2 = refine("33", "23")
    numpy.testing.assert_allclose(value_2, value, rtol=1e-4)
    numpy.testing.assert_allclose(sigma_2, sigma, rtol=1e-3)
