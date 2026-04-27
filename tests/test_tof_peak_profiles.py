import numpy

from cryspy.A_functions_base.powder_diffraction_tof import (
    calc_peak_shape_function,
    tof_non_convoluted_pseudo_voigt,
)
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile


def test_non_convoluted_pseudo_voigt_with_zero_gamma_is_gaussian():
    time = numpy.linspace(-3., 3., 41)
    time_hkl = numpy.array([0.25])
    sigma = numpy.array([0.9])
    gamma = numpy.array([0.])

    actual = tof_non_convoluted_pseudo_voigt(sigma, gamma, time, time_hkl)
    delta_t = time[:, numpy.newaxis] - time_hkl[numpy.newaxis, :]
    expected = numpy.exp(-0.5*numpy.square(delta_t/sigma)) / (
        numpy.sqrt(2.*numpy.pi)*sigma)

    numpy.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_non_convoluted_pseudo_voigt_uses_sigma_gamma_without_alpha_beta():
    time = numpy.linspace(-3., 3., 41)
    time_hkl = numpy.array([0.25])
    d = numpy.array([1.])
    sigmas = numpy.array([0.81, 0., 0.])
    gammas = numpy.array([0.35, 0., 0.])

    actual = calc_peak_shape_function(
        None, None, sigmas, d, time, time_hkl,
        gammas=gammas, peak_shape="non-conv-pseudo-Voigt")
    expected = tof_non_convoluted_pseudo_voigt(
        numpy.array([0.9]), numpy.array([0.35]), time, time_hkl)

    numpy.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_tof_profile_dictionary_non_convoluted_pseudo_voigt_has_no_alpha_beta():
    profile = TOFProfile(
        peak_shape="non-conv-pseudo-Voigt",
        sigma0=0.81,
        sigma1=0.,
        sigma2=0.,
        gamma0=0.35,
        gamma1=0.,
        gamma2=0.,
    )

    profile_dict = profile.get_dictionary()

    assert "profile_sigmas" in profile_dict
    assert "profile_gammas" in profile_dict
    assert "profile_alphas" not in profile_dict
    assert "profile_betas" not in profile_dict
