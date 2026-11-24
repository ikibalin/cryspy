"""Analytical derivatives of f_m_perp and of the polarized powder intensity
against finite differences (the keys once had typos: "eq_css", "flippe")."""
import numpy

from cryspy.A_functions_base.structure_factor import calc_f_m_perp_by_sft
from cryspy.A_functions_base.matrix_operations import calc_vector_product_v1_v2_v1
from cryspy.A_functions_base.integrated_intensity_powder_diffraction import \
    calc_powder_iint_1d_mix

STEP = 1e-6


def take_random_vectors(seed: int = 0, n_hkl: int = 4):
    rng = numpy.random.default_rng(seed)
    sft_ccs = rng.normal(size=(9, n_hkl)) + 1j * rng.normal(size=(9, n_hkl))
    eq_ccs = rng.normal(size=(3, n_hkl))
    eq_ccs = eq_ccs / numpy.linalg.norm(eq_ccs, axis=0)
    return sft_ccs, eq_ccs


def test_vector_product_v1_v2_v1_derivative_by_v1():
    rng = numpy.random.default_rng(1)
    v1 = rng.normal(size=(3, 5))
    v2 = rng.normal(size=(3, 5)) + 1j * rng.normal(size=(3, 5))
    res, dder = calc_vector_product_v1_v2_v1(v1, v2, flag_v1=True)
    assert dder["v1"].shape == (3, 3, 5)
    for j in range(3):
        v1_shifted = numpy.copy(v1)
        v1_shifted[j] += STEP
        numerical = (calc_vector_product_v1_v2_v1(v1_shifted, v2)[0] - res) / STEP
        assert numpy.allclose(dder["v1"][:, j], numerical, atol=1e-5)


def test_f_m_perp_by_sft_derivative_by_eq_ccs():
    sft_ccs, eq_ccs = take_random_vectors()
    magnetic_field = numpy.array([0., 0.3, 0.9], dtype=float)
    f_m_perp, dder = calc_f_m_perp_by_sft(
        sft_ccs, magnetic_field, eq_ccs, flag_eq_ccs=True)
    for j in range(3):
        eq_shifted = numpy.copy(eq_ccs)
        eq_shifted[j] += STEP
        numerical = (calc_f_m_perp_by_sft(sft_ccs, magnetic_field, eq_shifted)[0] - f_m_perp) / STEP
        assert numpy.allclose(dder["eq_ccs"][:, j], numerical, atol=1e-5)


def test_powder_iint_1d_mix_derivative_by_flipper():
    rng = numpy.random.default_rng(2)
    n_hkl = 3
    f_nucl = rng.normal(size=n_hkl) + 1j * rng.normal(size=n_hkl)
    f_m_perp = rng.normal(size=(3, n_hkl)) + 1j * rng.normal(size=(3, n_hkl))
    tensor_sigma = rng.normal(size=(9, n_hkl)) + 1j * rng.normal(size=(9, n_hkl))
    polarization, flipper, magnetic_field = 0.9, 0.95, 1.5

    def calc(flipper_efficiency, flag_flipper=False):
        return calc_powder_iint_1d_mix(
            f_nucl, f_m_perp, tensor_sigma, polarization, flipper_efficiency,
            magnetic_field, flag_flipper=flag_flipper)

    iint_plus, iint_minus, dder_plus, dder_minus = calc(flipper, flag_flipper=True)
    iint_plus_2, iint_minus_2 = calc(flipper + STEP)[:2]
    assert numpy.allclose(iint_plus_2, iint_plus)
    assert numpy.allclose(dder_minus["flipper"], (iint_minus_2 - iint_minus) / STEP, atol=1e-4)
