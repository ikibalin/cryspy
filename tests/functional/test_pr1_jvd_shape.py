"""Tests for the TOF Jorgensen-Von Dreele single-FWHM fix."""

import numpy as np

from cryspy.A_functions_base.powder_diffraction_tof import (
    tof_Jorgensen,
    tof_Jorgensen_VonDreele,
)


_NPTS = 401
_TIME = np.linspace(-200.0, 200.0, _NPTS)
_TIME_HKL = np.array([0.0])
_ALPHA = np.full(_NPTS, 0.5)
_BETA = np.full(_NPTS, 0.05)
_SIGMA = np.full(_NPTS, 20.0)


def test_jvd_reduces_to_gaussian_at_zero_gamma():
    gauss = tof_Jorgensen(_ALPHA, _BETA, _SIGMA, _TIME, _TIME_HKL)
    jvd = tof_Jorgensen_VonDreele(
        _ALPHA, _BETA, _SIGMA, np.zeros(_NPTS), _TIME, _TIME_HKL
    )
    np.testing.assert_allclose(jvd, gauss, rtol=1e-10, atol=1e-12)


def test_jvd_has_lorentzian_wings_for_nonzero_gamma():
    gauss = tof_Jorgensen(_ALPHA, _BETA, _SIGMA, _TIME, _TIME_HKL)
    jvd = tof_Jorgensen_VonDreele(
        _ALPHA, _BETA, _SIGMA, np.full(_NPTS, 15.0), _TIME, _TIME_HKL
    )
    far = np.abs(_TIME) > 120.0
    assert jvd[far].max() > 2.0 * gauss[far].max()


def test_jvd_golden_values_single_tch_fwhm():
    jvd = tof_Jorgensen_VonDreele(
        _ALPHA, _BETA, _SIGMA, np.full(_NPTS, 15.0), _TIME, _TIME_HKL
    )
    idx = [180, 200, 220, 260]
    expected = [0.00587542, 0.01106804, 0.01189444, 0.00408059]
    np.testing.assert_allclose(
        [jvd[i, 0] for i in idx], expected, rtol=1e-5
    )
