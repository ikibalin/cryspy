"""Extra broadening of magnetic reflections: Lorentzian width, hkl dependence, symmetry."""
import numpy
import pytest

from cryspy.A_functions_base.magnetic_broadening import (
    calc_aniso_by_isotropic_width,
    calc_aniso_from_independent,
    calc_independent_aniso_components,
    calc_principal_widths,
    calc_width_hkl,
    calc_y_mag_by_width,
    symmetrize_aniso,
)
from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_profile_pseudo_voight,
    calc_profile_pseudo_voight_extra_lorentz,
)

RAD = numpy.pi / 180.


def form_rotations(generators):
    """Closure of the group of integer 3x3 matrices with the inversion
    (coordinates x' = R x, acting on hkl as h' = R^T h)."""
    rotations = [numpy.eye(3, dtype=int), -numpy.eye(3, dtype=int)]
    flag_new = True
    while flag_new:
        flag_new = False
        for r_1 in list(rotations):
            for r_2 in list(generators):
                r_12 = numpy.matmul(r_1, numpy.asarray(r_2, dtype=int))
                if not any(numpy.array_equal(r_12, r) for r in rotations):
                    rotations.append(r_12)
                    flag_new = True
    return numpy.array(rotations, dtype=float)


C_4Z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
C_2X = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
C_2Y = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
C_3Z_HEX = [[0, -1, 0], [1, -1, 0], [0, 0, 1]]
C_3_111 = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]


def test_y_mag_is_the_fullprof_size_term():
    """FullProf: apparent size D = 360 lambda / (pi^2 Y) for the Lorentzian
    size term Y/cos(theta) (degrees); a half width w in Q is D = 2/w."""
    wavelength, width = 2.52, 0.02
    y_mag = calc_y_mag_by_width(width, wavelength)
    size_fullprof = 360. * wavelength / (numpy.pi ** 2 * y_mag)
    assert size_fullprof == pytest.approx(2. / width, rel=1e-12)
    assert calc_y_mag_by_width(0., wavelength) == 0.
    assert calc_y_mag_by_width(-width, wavelength) == pytest.approx(-y_mag)


def test_extra_lorentz_profile_is_the_profile_with_other_y():
    ttheta = numpy.linspace(10., 80., 3001) * RAD
    ttheta_hkl = numpy.array([20., 41.3, 70.]) * RAD
    args = (0.1, -0.08, 0.04, 0.002, 0.01, 0.03)
    u, v, w, i_g, x, y = args
    asymmetry = (0.05, 0.01, 0.02, 0.)
    for cutoff in (0., 8.):
        profile = calc_profile_pseudo_voight(
            ttheta, ttheta_hkl, *args, *asymmetry, cutoff_fwhm=cutoff)[0]
        profile_0 = calc_profile_pseudo_voight_extra_lorentz(
            ttheta, ttheta_hkl, *args, numpy.zeros(3), *asymmetry, cutoff_fwhm=cutoff)[0]
        assert numpy.allclose(profile_0, profile, rtol=0., atol=1e-13)
        for y_extra in (0.2, -0.02):
            profile_y = calc_profile_pseudo_voight(
                ttheta, ttheta_hkl, u, v, w, i_g, x, y + y_extra, *asymmetry,
                cutoff_fwhm=cutoff)[0]
            profile_extra = calc_profile_pseudo_voight_extra_lorentz(
                ttheta, ttheta_hkl, *args, numpy.full(3, y_extra), *asymmetry,
                cutoff_fwhm=cutoff)[0]
            assert numpy.allclose(profile_extra, profile_y, rtol=0., atol=1e-13)
    # different widths for the reflections
    profile_mixed = calc_profile_pseudo_voight_extra_lorentz(
        ttheta, ttheta_hkl, *args, numpy.array([0., 0.2, 0.]), *asymmetry)[0]
    profile_0 = calc_profile_pseudo_voight(ttheta, ttheta_hkl, *args, *asymmetry)[0]
    profile_y = calc_profile_pseudo_voight(
        ttheta, ttheta_hkl, u, v, w, i_g, x, y + 0.2, *asymmetry)[0]
    assert numpy.allclose(profile_mixed[:, 0], profile_0[:, 0])
    assert numpy.allclose(profile_mixed[:, 1], profile_y[:, 1])


def test_negative_width_removes_at_most_the_lorentzian_part():
    ttheta = numpy.linspace(30., 50., 2001) * RAD
    ttheta_hkl = numpy.array([40.]) * RAD
    args = (0., 0., 0.01, 0., 0., 0.05)
    gauss = calc_profile_pseudo_voight(ttheta, ttheta_hkl, 0., 0., 0.01, 0., 0., 0.,
                                       0., 0., 0., 0.)[0]
    profile = calc_profile_pseudo_voight_extra_lorentz(
        ttheta, ttheta_hkl, *args, numpy.array([-1.]), 0., 0., 0., 0.)[0]
    assert numpy.allclose(profile, gauss)


def test_half_width_in_q_is_the_width():
    """Independent check of the conversion: with a negligible resolution the
    magnetic peak is a Lorentzian in Q with the half width `width`."""
    wavelength, width = 2.52, 0.025
    ttheta_hkl = numpy.array([60. * RAD])
    ttheta = numpy.linspace(50., 70., 400001) * RAD
    y_mag = calc_y_mag_by_width(width, wavelength)
    profile = calc_profile_pseudo_voight_extra_lorentz(
        ttheta, ttheta_hkl, 0., 0., 1e-10, 0., 0., 0., numpy.array([y_mag]),
        0., 0., 0., 0.)[0][:, 0]
    q = 4. * numpy.pi * numpy.sin(0.5 * ttheta) / wavelength
    above = q[profile >= 0.5 * profile.max()]
    hwhm_q = 0.5 * (above.max() - above.min())
    # the width is taken at the point theta, so the peak is a Lorentzian in Q
    # up to the variation of cos(theta) over the peak (below 1 %)
    assert hwhm_q == pytest.approx(width, rel=1e-2)


@pytest.mark.parametrize("ucp_deg", [
    (5., 6., 7., 90., 90., 90.), (5., 6., 7., 90., 110., 90.),
    (3., 3., 5., 90., 90., 120.), (5., 6., 7., 80., 85., 100.)])
def test_isotropic_tensor_gives_the_same_width_for_every_direction(ucp_deg):
    ucp = numpy.array(ucp_deg) * numpy.array([1., 1., 1., RAD, RAD, RAD])
    for width in (0.03, -0.01):
        aniso = calc_aniso_by_isotropic_width(width, ucp)
        index_hkl = numpy.array([[1, 0, 0, 1, 1, 2, -1], [0, 1, 0, 1, -1, 1, 3],
                                 [0, 0, 1, 1, 2, -3, 1]], dtype=float)
        assert numpy.allclose(calc_width_hkl(index_hkl, ucp, aniso), width)
        assert numpy.allclose(calc_principal_widths(aniso, ucp)[0], width)


def test_width_along_reciprocal_axes():
    ucp = numpy.array([5., 6., 7., 90. * RAD, 90. * RAD, 90. * RAD])
    aniso = numpy.array([0.01, 0.02, 0.05, 0., 0., 0.])
    index_hkl = numpy.array([[2, 0, 0], [0, 1, 0], [0, 0, 3]], dtype=float)
    assert numpy.allclose(calc_width_hkl(index_hkl, ucp, aniso), [0.01, 0.02, 0.05])
    # (1 1 0) is in between, weighted by the squared direction cosines
    width_110 = calc_width_hkl(numpy.array([[1.], [1.], [0.]]), ucp, aniso)[0]
    cos_sq_a = (1. / 5.) ** 2 / ((1. / 5.) ** 2 + (1. / 6.) ** 2)
    assert width_110 == pytest.approx(0.01 * cos_sq_a + 0.02 * (1. - cos_sq_a))
    widths, axes = calc_principal_widths(aniso, ucp)
    assert numpy.allclose(widths, [0.01, 0.02, 0.05])
    assert numpy.allclose(numpy.abs(axes[:, 2]), [0., 0., 1.])


@pytest.mark.parametrize("generators, ucp_deg, expected_independent", [
    ([C_4Z, C_2X], (4., 4., 9., 90., 90., 90.), [True, False, True, False, False, False]),
    ([C_3Z_HEX], (4., 4., 9., 90., 90., 120.), [True, False, True, False, False, False]),
    ([C_2X, C_2Y], (4., 5., 9., 90., 90., 90.), [True, True, True, False, False, False]),
    ([C_2Y], (4., 5., 9., 90., 105., 90.), [True, True, True, False, True, False]),
    ([C_3_111, C_4Z], (4., 4., 4., 90., 90., 90.), [True, False, False, False, False, False]),
    ([], (4., 5., 9., 80., 85., 100.), [True] * 6),
])
def test_symmetry_of_the_tensor(generators, ucp_deg, expected_independent):
    ucp = numpy.array(ucp_deg) * numpy.array([1., 1., 1., RAD, RAD, RAD])
    rotations = form_rotations(generators)
    projection, independent = calc_independent_aniso_components(rotations, ucp)
    assert list(independent) == expected_independent
    numpy.random.seed(5)
    aniso = numpy.random.uniform(-0.01, 0.04, 6)
    aniso_sym = symmetrize_aniso(aniso, rotations, ucp)
    assert numpy.allclose(projection @ aniso, aniso_sym)
    # symmetric tensors keep their values; the dependent components are
    # restored from the independent ones
    assert numpy.allclose(symmetrize_aniso(aniso_sym, rotations, ucp), aniso_sym)
    trial = numpy.where(independent, aniso_sym, 0.)
    assert numpy.allclose(calc_aniso_from_independent(trial, projection, independent), aniso_sym)
    # symmetry-equivalent reflections have the same width
    index_hkl = numpy.array([[1.], [2.], [3.]])
    width_0 = calc_width_hkl(index_hkl, ucp, aniso_sym)
    for rotation in rotations:
        index_eq = rotation.transpose() @ index_hkl
        assert numpy.allclose(calc_width_hkl(index_eq, ucp, aniso_sym), width_0)


def test_hexagonal_isotropic_tensor_has_the_u_ij_relations():
    ucp = numpy.array([4., 4., 9., 90. * RAD, 90. * RAD, 120. * RAD])
    aniso = calc_aniso_by_isotropic_width(0.02, ucp)
    # U_12 = U_11 / 2 for the hexagonal axes (angle a*, b* = 60 degrees)
    assert aniso[3] == pytest.approx(0.5 * aniso[0])
    rotations = form_rotations([C_3Z_HEX])
    assert numpy.allclose(symmetrize_aniso(aniso, rotations, ucp), aniso)
