import numpy
import scipy.special

from cryspy.A_functions_base.powder_diffraction_tof import (
    calc_d_min_max_by_time_thermal_neutrons,
    calc_hpv_eta,
    tof_Jorgensen,
    tof_Jorgensen_VonDreele,
)
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters


def test_calc_d_min_max_by_time_thermal_neutrons_supports_zero_dtt2():
    time = numpy.array([3000.0, 19000.0], dtype=float)
    zero = numpy.array([2.921], dtype=float)
    dtt1 = numpy.array([6167.247], dtype=float)
    dtt2 = numpy.array([0.0], dtype=float)

    d_min_max = calc_d_min_max_by_time_thermal_neutrons(time, zero, dtt1, dtt2)

    expected = numpy.array([
        (time.min() - zero[0]) / dtt1[0],
        (time.max() - zero[0]) / dtt1[0],
    ], dtype=float)
    numpy.testing.assert_allclose(d_min_max, expected)

    tof_parameters = TOFParameters(
        zero=float(zero[0]),
        dtt1=float(dtt1[0]),
        dtt2=float(dtt2[0]),
        ttheta_bank=145.0,
    )
    numpy.testing.assert_allclose(tof_parameters.calc_d_min_max(time), expected)


def _expected_lorentz_profile(alpha, beta, gamma, time, time_hkl):
    norm = 0.5 * alpha * beta / (alpha + beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d - time_hkl_2d

    z1_2d = alpha[:, numpy.newaxis] * delta_2d + (0.5j * alpha * gamma)[:, numpy.newaxis]
    z2_2d = -beta[:, numpy.newaxis] * delta_2d + (0.5j * beta * gamma)[:, numpy.newaxis]

    with numpy.errstate(over="ignore", invalid="ignore"):
        term_1 = numpy.exp(z1_2d) * scipy.special.exp1(z1_2d)
        term_2 = numpy.exp(z2_2d) * scipy.special.exp1(z2_2d)

    term_1[numpy.isnan(term_1)] = 0.0
    term_1[numpy.isinf(term_1)] = 0.0
    term_2[numpy.isnan(term_2)] = 0.0
    term_2[numpy.isinf(term_2)] = 0.0

    return -2.0 * norm[:, numpy.newaxis] * (
        numpy.imag(term_1) + numpy.imag(term_2)
    ) / numpy.pi


def test_tof_jorgensen_von_dreele_matches_expected_lorentz_term():
    alpha = numpy.array([0.4], dtype=float)
    beta = numpy.array([0.2], dtype=float)
    sigma = numpy.array([10.0], dtype=float)
    gamma = numpy.array([5.0], dtype=float)
    time = numpy.linspace(-50.0, 50.0, 11)
    time_hkl = numpy.array([0.0], dtype=float)

    actual = tof_Jorgensen_VonDreele(alpha, beta, sigma, gamma, time, time_hkl)

    gaussian = tof_Jorgensen(alpha, beta, sigma, time, time_hkl)
    _, eta = calc_hpv_eta(sigma, gamma)
    lorentz = _expected_lorentz_profile(alpha, beta, gamma, time, time_hkl)
    expected = (1.0 - eta)[:, numpy.newaxis] * gaussian + eta[:, numpy.newaxis] * lorentz

    numpy.testing.assert_allclose(actual, expected)
    assert numpy.isfinite(actual).all()
    assert numpy.all(actual >= 0.0)
