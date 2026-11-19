"""Texture by spherical harmonics: comparison with GSAS-II and invariants."""
import numpy
import pytest

from cryspy.A_functions_base.absorption import calc_rotation_matrix_busing_levy
from cryspy.A_functions_base.texture_spherical_harmonics import (
    LAUE_CLASSES,
    calc_crystal_angles,
    calc_crystal_harmonics,
    calc_laue_class_by_rotations,
    calc_sample_angles,
    calc_texture_index,
    calc_texture_sh_correction,
    form_texture_sh_lmn,
)

RAD = numpy.pi / 180.

# GSAS-II 2DTexture tutorial (NiTi-A.gpx), phase B2 (Pm-3m, a = 3.01 Å),
# texture of the phase, order 12, cylindrical; Debye-Scherrer, omega = 0,
# chi = 90, phi = 0. Columns: h, k, l, 2theta (deg), texture correction
# of GSAS-II (Prfo, stored as float32).
NITI_B2_COEFFICIENTS = {(4, 0, 1): -3.199470401782551, (6, 0, 1): -3.4035880044915556,
                        (8, 0, 1): 1.6503831440459107, (10, 0, 1): 2.109431163504195,
                        (12, 0, 1): -0.5259836645205405, (12, 0, 2): 0.848140828633415}
NITI_B2_PRFO = {
    0.: [(1, 0, 0, 3.37798, 0.2120605707168579), (1, 1, 0, 4.77788, 2.150869131088257),
         (1, 1, 1, 5.85253, 0.788860023021698), (2, 0, 0, 6.7589, 0.2120605707168579),
         (2, 1, 0, 7.55778, 0.7141231298446655), (2, 1, 1, 8.28033, 1.1046640872955322)],
    10.: [(1, 0, 0, 3.37369, 0.18190287053585052), (1, 1, 0, 4.77181, 1.360916256904602),
          (1, 1, 1, 5.84509, 1.195592999458313), (2, 0, 0, 6.75031, 0.18178175389766693),
          (2, 1, 0, 7.54817, 1.070197582244873), (2, 1, 1, 8.2698, 1.0095504522323608)],
}

UCP_CUBIC = numpy.array([3.01, 3.01, 3.01, 90. * RAD, 90. * RAD, 90. * RAD])


@pytest.mark.parametrize("azimuth", sorted(NITI_B2_PRFO))
def test_correction_as_gsas2_niti(azimuth):
    table = numpy.array(NITI_B2_PRFO[azimuth])
    index_hkl = table[:, :3].transpose()
    lmn = numpy.array(list(NITI_B2_COEFFICIENTS)).transpose()
    coefficient = numpy.array(list(NITI_B2_COEFFICIENTS.values()))
    psi, gamma = calc_sample_angles(0.5 * table[:, 3] * RAD, 0., 90. * RAD, 0., azimuth * RAD)
    correction = calc_texture_sh_correction(index_hkl, UCP_CUBIC, "m-3m", "", "cylindrical",
                                            lmn, coefficient, psi, gamma)[0]
    # GSAS-II stores the correction as float32
    assert numpy.allclose(correction, table[:, 4], rtol=1e-6, atol=1e-6)


def test_number_of_coefficients_as_gsas2():
    # numbers of terms of GSAS-II GenSHCoeffT for the tutorial phases
    assert form_texture_sh_lmn("m-3m", "cylindrical", 12).shape == (3, 6)
    assert form_texture_sh_lmn("2/m", "cylindrical", 10).shape == (3, 35)
    lmn = form_texture_sh_lmn("mmm", "mmm", 4)
    assert set(map(tuple, lmn.transpose())) == {
        (l, m, n) for l in (2, 4) for m in range(0, l + 1, 2) for n in range(0, l + 1, 2)}


def test_crystal_angles_as_gsas2():
    # GSAS-II CrsAng for orthogonal axes: polar angle from c*, azimuth from a*
    ucp = numpy.array([4.1, 5.3, 7.7, 90. * RAD, 90. * RAD, 90. * RAD])
    index_hkl = numpy.array([[1, 2, 3], [2, -1, 1], [0, 1, 0]], dtype=float).transpose()
    polar, azimuth = calc_crystal_angles(index_hkl, ucp, "mmm")
    q = index_hkl / ucp[:3, None]
    assert numpy.allclose(numpy.cos(polar), q[2] / numpy.linalg.norm(q, axis=0))
    assert numpy.allclose(azimuth, numpy.arctan2(q[1], q[0]))
    # hexagonal axes: x || a, a* is 30 degrees and a* + b* 60 degrees from it
    ucp_hex = numpy.array([3., 3., 5., 90. * RAD, 90. * RAD, 120. * RAD])
    polar, azimuth = calc_crystal_angles([[1, 1], [0, 1], [0, 0]], ucp_hex, "6/mmm")
    assert numpy.allclose(polar, 0.5 * numpy.pi)
    assert numpy.allclose(numpy.degrees(azimuth), [30., 60.])
    # rhombohedral axes: (1 1 1) along z
    ucp_rh = numpy.array([5., 5., 5., 70. * RAD, 70. * RAD, 70. * RAD])
    polar = calc_crystal_angles([[1], [1], [1]], ucp_rh, "-3mR")[0]
    assert numpy.allclose(polar, 0.)


def form_rotations(generators):
    """Closure of the group of integer 3x3 matrices with the inversion."""
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
    return numpy.array(rotations)


C_4Z = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
C_2X = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
C_2Y = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
C_3Z_HEX = [[0, -1, 0], [1, -1, 0], [0, 0, 1]]
C_2Z = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
C_2_110 = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
C_2_1M10 = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
C_3_111 = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
C_6Z_HEX = [[1, -1, 0], [1, 0, 0], [0, 0, 1]]


@pytest.mark.parametrize("generators, ucp_deg, expected", [
    ([], (3., 4., 5., 80., 85., 95.), ("-1", "")),
    ([C_2Y], (3., 4., 5., 90., 100., 90.), ("2/m", "b")),
    ([C_2Z], (3., 4., 5., 90., 90., 100.), ("2/m", "c")),
    ([C_2X, C_2Y], (3., 4., 5., 90., 90., 90.), ("mmm", "")),
    ([C_4Z], (3., 3., 5., 90., 90., 90.), ("4/m", "")),
    ([C_4Z, C_2X], (3., 3., 5., 90., 90., 90.), ("4/mmm", "")),
    ([C_3Z_HEX], (3., 3., 5., 90., 90., 120.), ("-3", "")),
    ([C_3Z_HEX, C_2_110], (3., 3., 5., 90., 90., 120.), ("-3m1", "")),
    ([C_3Z_HEX, C_2_1M10], (3., 3., 5., 90., 90., 120.), ("-31m", "")),
    ([C_3_111, C_2_1M10], (5., 5., 5., 70., 70., 70.), ("-3mR", "")),
    ([C_6Z_HEX], (3., 3., 5., 90., 90., 120.), ("6/m", "")),
    ([C_6Z_HEX, C_2_110], (3., 3., 5., 90., 90., 120.), ("6/mmm", "")),
    ([C_3_111], (5., 5., 5., 70., 70., 70.), ("-3R", "")),
    ([C_3_111, C_4Z], (4., 4., 4., 90., 90., 90.), ("m-3m", "")),
    ([C_3_111, C_2X], (4., 4., 4., 90., 90., 90.), ("m-3", "")),
])
def test_laue_class_by_rotations(generators, ucp_deg, expected):
    # matrices of the coordinates x' = R x, acting on hkl as h' = R^T h
    rotations = form_rotations(generators)
    ucp = numpy.array(ucp_deg, dtype=float)
    ucp[3:] *= RAD
    assert calc_laue_class_by_rotations(rotations, ucp) == expected


@pytest.mark.parametrize("laue_class", [h for h in LAUE_CLASSES if h not in ("-3R", "-3mR")])
def test_crystal_harmonics_are_invariant(laue_class):
    # the harmonics are the same for the reflections equivalent by the 4-,
    # 3- or 6-fold axis along c
    ucp_deg = (3., 3., 5., 90., 90., 120.) if laue_class[:2] in ("-3", "6/") \
        else (4., 4., 4., 90., 90., 90.)
    ucp = numpy.array(ucp_deg) * numpy.array([1, 1, 1, RAD, RAD, RAD])
    rotation = {"4": C_4Z, "-": C_3Z_HEX, "6": C_6Z_HEX, "m": C_3_111}.get(laue_class[0])
    if laue_class in ("-1", "2/m", "mmm"):
        rotation = -numpy.eye(3)
    lmn = form_texture_sh_lmn(laue_class, "cylindrical", 8)
    index_hkl = numpy.random.default_rng(3).uniform(-2., 2., (3, 7))
    index_hkl_rot = numpy.matmul(numpy.transpose(rotation), index_hkl)
    k_1 = calc_crystal_harmonics(lmn[[0, 2]], laue_class,
                                 *calc_crystal_angles(index_hkl, ucp, laue_class, "b"))
    k_2 = calc_crystal_harmonics(lmn[[0, 2]], laue_class,
                                 *calc_crystal_angles(index_hkl_rot, ucp, laue_class, "b"))
    assert numpy.allclose(k_1, k_2, atol=1e-10)


def test_sample_angles_are_busing_levy():
    # q in the laboratory system of Busing & Levy (y along the beam, z up,
    # detector azimuth from +x towards +z), q_phi = R^T q_lab, R = Omega Chi Phi;
    # the texture axis is x of the phi frame, the azimuth from -z towards +y
    rng = numpy.random.default_rng(11)
    theta = rng.uniform(1., 80., 25) * RAD
    azimuth = rng.uniform(-180., 180., 25) * RAD
    q_lab = numpy.stack([numpy.cos(theta) * numpy.cos(azimuth), -numpy.sin(theta),
                         numpy.cos(theta) * numpy.sin(azimuth)])
    for omega, chi, phi in [(0., 0., 0.), (0., 90., 0.), (23., 41., -67.), (-110., 15., 200.)]:
        matrix_r = calc_rotation_matrix_busing_levy([omega], [chi], [phi])[:, :, 0]
        q_phi = numpy.matmul(matrix_r.transpose(), q_lab)
        psi, gamma = calc_sample_angles(theta, omega * RAD, chi * RAD, phi * RAD, azimuth)
        assert numpy.allclose(numpy.cos(psi), q_phi[0])
        assert numpy.allclose(gamma, numpy.arctan2(q_phi[1], -q_phi[2]))


def test_sample_angles_default_axis_is_vertical():
    # chi = 90: the texture axis is vertical; q of the detector azimuth 90
    # (vertical scattering plane) makes the angle theta with the vertical
    theta = numpy.array([10., 30.]) * RAD
    psi = calc_sample_angles(theta, 0., 90. * RAD, 0., 0.5 * numpy.pi)[0]
    assert numpy.allclose(numpy.abs(numpy.cos(psi)), numpy.cos(theta))
    psi = calc_sample_angles(theta, 0., 90. * RAD, 0., 0.)[0]
    assert numpy.allclose(psi, 0.5 * numpy.pi)


def test_random_powder_and_derivatives():
    lmn = form_texture_sh_lmn("4/mmm", "mmm", 6)
    ucp = numpy.array([3., 3., 5., 90. * RAD, 90. * RAD, 90. * RAD])
    rng = numpy.random.default_rng(5)
    index_hkl = rng.integers(-3, 4, (3, 9))
    index_hkl[2, 0] = 1
    psi, gamma = rng.uniform(0., numpy.pi, (4, 9)), rng.uniform(-numpy.pi, numpy.pi, (4, 9))
    correction, dder = calc_texture_sh_correction(
        index_hkl, ucp, "4/mmm", "", "mmm", lmn, numpy.zeros(lmn.shape[1]), psi, gamma,
        flag_texture_sh_coefficient=True)
    assert numpy.allclose(correction, 1.)
    coefficient = rng.normal(size=lmn.shape[1])
    correction_2 = calc_texture_sh_correction(
        index_hkl, ucp, "4/mmm", "", "mmm", lmn, coefficient, psi, gamma)[0]
    # the correction is linear in the coefficients
    der = dder["texture_sh_coefficient"]
    assert der.shape == psi.shape + (lmn.shape[1],)
    assert numpy.allclose(correction_2, 1. + numpy.matmul(der, coefficient))


def test_mean_over_sample_directions_is_one():
    # the integral of every harmonic of l > 0 over the sphere is zero, so the
    # texture does not change the mean intensity
    lmn = form_texture_sh_lmn("m-3m", "-1", 8)
    coefficient = numpy.random.default_rng(1).normal(size=lmn.shape[1])
    cos_psi, w_psi = numpy.polynomial.legendre.leggauss(30)
    gamma = numpy.arange(40) * 2. * numpy.pi / 40
    psi, gamma = numpy.meshgrid(numpy.arccos(cos_psi), gamma, indexing="ij")
    weight = w_psi[:, None] * numpy.full(psi.shape, 2. * numpy.pi / 40) / (4. * numpy.pi)
    index_hkl = numpy.array([[1, 2, 3], [0, 1, 4]], dtype=float).transpose()
    correction = calc_texture_sh_correction(
        index_hkl, UCP_CUBIC, "m-3m", "", "-1", lmn, coefficient,
        psi[..., None] + 0. * index_hkl[0], gamma[..., None] + 0. * index_hkl[0])[0]
    assert numpy.allclose(numpy.einsum("ijh,ij->h", correction, weight), 1., atol=1e-10)


def test_texture_index():
    lmn = numpy.array([[2, 4], [0, 0], [0, 1]])
    assert calc_texture_index(lmn, [0., 0.]) == pytest.approx(1.)
    assert calc_texture_index(lmn, [5., 3.]) == pytest.approx(1. + 25. / 5. + 9. / 9.)


def test_unknown_symmetry():
    with pytest.raises(ValueError):
        form_texture_sh_lmn("mmm", "cubic", 4)
