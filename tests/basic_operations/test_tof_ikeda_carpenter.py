"""Ikeda-Carpenter pulse convoluted with a pseudo-Voigt (issue #39).

The closed form is checked against a direct numerical convolution of the
Ikeda-Carpenter function (Nucl. Instrum. Methods A 239, 536 (1985)) with
the pseudo-Voigt, so the reference does not use the code under test.
"""
import numpy
import pytest
import scipy.integrate

from cryspy.A_functions_base.powder_diffraction_tof import (
    calc_hpv_eta,
    calc_peak_shape_function,
    inv_8ln2,
)
from cryspy.A_functions_base.powder_diffraction_tof_ikeda_carpenter import (
    calc_ikeda_carpenter_alpha_beta_r,
    tof_ikeda_carpenter,
)
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile


def ikeda_carpenter(tau, alpha, beta, r):
    a_b = alpha - beta
    value = 0.5 * alpha * (
        (1. - r) * (alpha * tau)**2 * numpy.exp(-alpha * tau)
        + 2. * r * alpha**2 * beta / a_b**3 * (
            numpy.exp(-beta * tau) - numpy.exp(-alpha * tau) * (
                1. + a_b * tau + 0.5 * (a_b * tau)**2)))
    return numpy.where(tau >= 0., value, 0.)


def pseudo_voigt(delta, sigma, gamma):
    h_pv, eta = calc_hpv_eta(numpy.array([sigma * numpy.sqrt(8. * numpy.log(2.))]),
                             numpy.array([gamma]))
    h_pv, eta = h_pv[0], eta[0]
    sigma_c = h_pv * numpy.sqrt(inv_8ln2)
    gauss = numpy.exp(-0.5 * (delta / sigma_c)**2) / (sigma_c * numpy.sqrt(2. * numpy.pi))
    lorentz = 0.5 * h_pv / numpy.pi / (delta**2 + 0.25 * h_pv**2)
    return (1. - eta) * gauss + eta * lorentz


def convolute(delta, alpha, beta, r, sigma, gamma):
    def integrand(tau):
        return ikeda_carpenter(tau, alpha, beta, r) * pseudo_voigt(delta - tau, sigma, gamma)
    split = max(delta, 0.) + 1.
    return (scipy.integrate.quad(integrand, 0., split, limit=1000, epsabs=1e-14)[0] +
            scipy.integrate.quad(integrand, split, numpy.inf, limit=1000, epsabs=1e-14)[0])


@pytest.mark.parametrize("alpha, beta, r, sigma, gamma", [
    (0.5, 0.03, 0.4, 8., 3.),     # usual moderator, both components
    (0.2, 0.05, 0.0, 3., 0.),     # no slow term, Gaussian only
    (1.0, 0.02, 0.9, 20., 10.),   # mostly slow term, wide Lorentzian
    (0.05, 0.04, 0.5, 5., 2.),    # alpha close to beta
    (2.0, 0.03, 0.4, 100., 30.),  # alpha sigma >> 1 (cancellations)
])
def test_profile_equals_numerical_convolution(alpha, beta, r, sigma, gamma):
    delta = numpy.array([-800., -60., -10., -2., 0., 3., 10., 40., 150., 600., 1200.])
    n_point = delta.size
    profile = tof_ikeda_carpenter(
        numpy.full(n_point, alpha), numpy.full(n_point, beta), numpy.full(n_point, r),
        numpy.full(n_point, sigma), numpy.full(n_point, gamma),
        delta, numpy.array([0.]))[:, 0]
    reference = numpy.array([convolute(d, alpha, beta, r, sigma, gamma) for d in delta])
    # 1e-6 of the maximum: the quadrature and the erfc cancellations for
    # alpha sigma = 200 are at 1e-7
    numpy.testing.assert_allclose(profile, reference, rtol=0., atol=1e-6 * reference.max())


def test_profile_has_unit_area_and_tail_to_longer_time():
    time = numpy.linspace(-500., 6000., 130001)
    n_point = time.size
    profile = tof_ikeda_carpenter(
        numpy.full(n_point, 0.4), numpy.full(n_point, 0.03), numpy.full(n_point, 0.5),
        numpy.full(n_point, 5.), numpy.zeros(n_point), time, numpy.array([0.]))[:, 0]
    assert numpy.trapezoid(profile, time) == pytest.approx(1., abs=1e-6)
    # rounding of the cancelling terms only
    assert profile.min() > -1e-12 * profile.max()
    # the moderator delays neutrons: more weight after the peak position
    assert profile[time > 0.].sum() > 3. * profile[time < 0.].sum()


def test_far_points_use_asymptotic_lorentzian():
    # |alpha delta| > 600: exp(z) E1(z) by its asymptotic series
    time = numpy.array([-20000., -3000., 3000., 20000.])
    n_point = time.size
    profile = tof_ikeda_carpenter(
        numpy.full(n_point, 0.5), numpy.full(n_point, 0.03), numpy.full(n_point, 0.4),
        numpy.full(n_point, 4.), numpy.full(n_point, 3.), time, numpy.array([0.]))[:, 0]
    assert numpy.all(numpy.isfinite(profile))
    reference = numpy.array([convolute(d, 0.5, 0.03, 0.4, 4., 3.) for d in time])
    numpy.testing.assert_allclose(profile, reference, rtol=1e-4)


def test_alpha_beta_r_follow_mantid_parametrisation():
    wavelength = numpy.array([0.5, 1., 2.])
    alpha, beta, r = calc_ikeda_carpenter_alpha_beta_r(wavelength, 1.6, 1.5, 31.9, 46.)
    numpy.testing.assert_allclose(alpha, 1. / (1.6 + 1.5 * wavelength))
    numpy.testing.assert_allclose(beta, 1. / 31.9)
    numpy.testing.assert_allclose(r, numpy.exp(-81.799 / (46. * wavelength**2)))
    assert numpy.all(calc_ikeda_carpenter_alpha_beta_r(wavelength, 1.6, 1.5, 31.9, 0.)[2] == 0.)


def test_peak_shape_function_with_cutoff():
    time = numpy.linspace(4000., 6000., 801)
    time_hkl = numpy.array([4500., 5200.])
    d = time / 6167.
    kwargs = dict(gammas=numpy.array([0., 2., 0.]), peak_shape="Ikeda-Carpenter",
                  kappa=46., ttheta_bank=numpy.deg2rad(145.))
    args = (numpy.array([1.6, 1.5]), numpy.array([31.9]), numpy.array([10., 30., 0.]),
            d, time, time_hkl)
    full = calc_peak_shape_function(*args, **kwargs)
    cut = calc_peak_shape_function(*args, cutoff_fwhm=20., **kwargs)
    inside = cut != 0.
    numpy.testing.assert_allclose(cut[inside], full[inside])
    assert full.shape == (time.size, time_hkl.size)
    assert numpy.all(numpy.isfinite(full))


def test_tof_profile_item_ikeda_carpenter():
    profile = TOFProfile(
        peak_shape="Ikeda-Carpenter", sigma0=10., sigma1=30., sigma2=0.,
        gamma0=0., gamma1=2., gamma2=0., alpha0=1.6, alpha1=1.5, beta0=31.9, kappa=46.)
    profile.kappa_refinement = True
    dict_profile = profile.get_dictionary()
    numpy.testing.assert_allclose(dict_profile["profile_alphas"], [1.6, 1.5])
    numpy.testing.assert_allclose(dict_profile["profile_betas"], [31.9])
    numpy.testing.assert_allclose(dict_profile["profile_kappa"], [46.])
    assert dict_profile["flags_profile_kappa"].tolist() == [True]
    text = profile.to_cif()
    assert "_tof_profile_kappa" in text and "_tof_profile_beta1" not in text


def test_tof_profile_item_ikeda_carpenter_needs_kappa():
    profile = TOFProfile(
        peak_shape="Ikeda-Carpenter", sigma0=10., sigma1=30., sigma2=0.,
        gamma0=0., gamma1=2., gamma2=0., alpha0=1.6, alpha1=1.5, beta0=31.9)
    with pytest.raises(ValueError, match="kappa"):
        profile.get_dictionary()
