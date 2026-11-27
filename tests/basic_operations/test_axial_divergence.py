"""Finger-Cox-Jephcoat axial divergence: shifts and weights of the quadrature."""
import numpy
import pytest
from scipy.integrate import quad

from cryspy.A_functions_base.powder_diffraction_axial_divergence import (
    calc_axial_divergence_shifts, calc_profile_with_axial_divergence)


def calc_moments_fcj(ttheta: float, s_l: float, d_l: float):
    """Mean and rms shift by direct integration of the weight of FCJ (1994),
    eq. 4-7, over 2phi (an independent route: singular weight, scipy quad)."""
    cos_ttheta = numpy.cos(ttheta)

    def weight(ttheta_phi):
        h_l = numpy.sqrt(max(numpy.cos(ttheta_phi) ** 2 / cos_ttheta ** 2 - 1., 0.))
        if h_l <= 0. or h_l >= s_l + d_l:
            return 0.
        a = 2. * min(s_l, d_l) if h_l < abs(s_l - d_l) else s_l + d_l - h_l
        return a / (h_l * abs(numpy.cos(ttheta_phi)))

    ttheta_min = numpy.arccos(numpy.clip(cos_ttheta * numpy.sqrt(1. + (s_l + d_l) ** 2), -1., 1.))
    limits = sorted([ttheta_min, ttheta])
    norm = quad(weight, *limits, limit=500)[0]
    mean = quad(lambda p: (p - ttheta) * weight(p), *limits, limit=500)[0] / norm
    second = quad(lambda p: (p - ttheta) ** 2 * weight(p), *limits, limit=500)[0] / norm
    return mean, numpy.sqrt(second)


@pytest.mark.parametrize("ttheta_deg", [10., 30., 60., 120., 160.])
@pytest.mark.parametrize("s_l, d_l", [(0.02, 0.03), (0.04, 0.01), (0.025, 0.025)])
def test_shifts_against_direct_integration(ttheta_deg, s_l, d_l):
    ttheta = numpy.radians(ttheta_deg)
    shift, weight = calc_axial_divergence_shifts(numpy.array([ttheta]), s_l, d_l)
    mean, rms = calc_moments_fcj(ttheta, s_l, d_l)
    assert weight.sum() == pytest.approx(1.)
    assert (shift * weight).sum() == pytest.approx(mean, rel=1e-3)
    assert numpy.sqrt((shift ** 2 * weight).sum()) == pytest.approx(rms, rel=1e-3)


def test_sign_of_the_asymmetry_and_no_shift():
    shift, weight = calc_axial_divergence_shifts(numpy.radians([20., 90., 150.]), 0.02, 0.03)
    mean = (shift * weight).sum(axis=0)
    assert mean[0] < 0. and mean[2] > 0.
    assert mean[1] == pytest.approx(0., abs=1e-12)
    shift, weight = calc_axial_divergence_shifts(numpy.radians([20., 150.]), 0., 0.03)
    numpy.testing.assert_allclose(shift, 0.)
    numpy.testing.assert_allclose(weight, 1.)


def test_profile_keeps_area_and_ttheta_max():
    ttheta = numpy.radians(numpy.linspace(5., 70., 6501))
    ttheta_hkl = numpy.radians([15., 50.])
    width = numpy.radians(0.2)

    def gauss(ttheta_shifted, flags_hkl):
        z = (ttheta[:, None] - ttheta_shifted[None, :]) / width
        return numpy.exp(-0.5 * z ** 2) / (width * numpy.sqrt(2. * numpy.pi))

    profile = calc_profile_with_axial_divergence(
        gauss, ttheta_hkl, 0.02, 0.03, ttheta_max=numpy.radians(30.))
    step = ttheta[1] - ttheta[0]
    numpy.testing.assert_allclose(profile.sum(axis=0) * step, 1., rtol=1e-4)
    # the peak below 2theta_max moves to low angles, the other one stays
    centroid = (ttheta[:, None] * profile).sum(axis=0) / profile.sum(axis=0)
    assert centroid[0] < ttheta_hkl[0] - numpy.radians(0.01)
    numpy.testing.assert_allclose(profile[:, 1], gauss(ttheta_hkl, None)[:, 1])
