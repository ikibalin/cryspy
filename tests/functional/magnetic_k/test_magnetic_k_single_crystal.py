"""Unpolarized single-crystal intensities of a structure with a propagation
vector, given by fractional indices H = h ± k, are compared with the same
structure written explicitly in the magnetic supercell (independent mcif code
for ordered moments, integer indices in the supercell).

Magnetic domains are the twins of ``jana_twin``: the intensity of the
k-structure with two domains must be the fraction-weighted sum of two
supercell calculations at H and at M·H.
"""
import numpy
import pytest

import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_calc_chi_sq_by_dictionary,
)
import test_magnetic_k_supercell as supercell

DOMAIN_MATRIX = numpy.array([[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
DOMAIN_FRACTION = 0.3

DIFFRN = """
data_sc

_setup_wavelength 1.2
_setup_field 0.0

_diffrn_radiation_polarization 0.0
_diffrn_radiation_efficiency 1.0

_diffrn_orient_matrix_ub_11 0.2
_diffrn_orient_matrix_ub_12 0.0
_diffrn_orient_matrix_ub_13 0.0
_diffrn_orient_matrix_ub_21 0.0
_diffrn_orient_matrix_ub_22 0.2
_diffrn_orient_matrix_ub_23 0.0
_diffrn_orient_matrix_ub_31 0.0
_diffrn_orient_matrix_ub_32 0.0
_diffrn_orient_matrix_ub_33 0.2

_phase_label phase
_phase_scale 1.0
{twin}
loop_
_diffrn_refln_index_h
_diffrn_refln_index_k
_diffrn_refln_index_l
_diffrn_refln_intensity
_diffrn_refln_intensity_sigma
{refln}
"""


def strip_experiment(rcif: str) -> str:
    return rcif[:rcif.index("data_pd")]


def form_diffrn(index_hkl, twin_matrix=None, fraction=None) -> str:
    refln = "\n".join(f"{h:.8f} {k:.8f} {l:.8f} 1.0 1.0"
                      for h, k, l in index_hkl.transpose())
    twin = ""
    if twin_matrix is not None:
        twin = ("\nloop_\n_jana_twin_id\n" + "".join(
            f"_jana_twin_matrix[{i}][{j}]\n" for i in (1, 2, 3) for j in (1, 2, 3))
            + "_jana_twin_fraction\n2 " + " ".join(
                f"{v:.8f}" for v in twin_matrix.flatten()) + f" {fraction}\n")
    return DIFFRN.format(twin=twin, refln=refln)


def calc_intensity(rcif: str):
    obj = cryspy.str_to_globaln(rcif)
    dict_obj = obj.get_dictionary()
    dict_in_out = {}
    rhochi_calc_chi_sq_by_dictionary(
        dict_obj, dict_in_out=dict_in_out, flag_use_precalculated_data=False,
        flag_calc_analytical_derivatives=False)
    return dict_in_out["diffrn_sc"]["intensity_calc"]


def form_index_hkl(k_vector):
    """Nuclear reflections, satellites of +k and of -k, and one index that
    is neither (zero intensity)."""
    index_h = numpy.array([[0, 0, 1], [1, 1, 0], [0, 2, 0], [1, 0, 1],
                           [2, 1, 1], [0, 1, 2], [1, 2, -1]], dtype=float).transpose()
    k_vector = numpy.asarray(k_vector)[:, numpy.newaxis]
    return numpy.concatenate(
        [index_h, index_h + k_vector, index_h - k_vector,
         numpy.array([[0.25], [0.], [0.]])], axis=1)


def check_with_supercell(case_name: str):
    case = supercell.CASES[case_name]
    index_hkl = form_index_hkl(case["k_vector"])
    rcif_k = strip_experiment(supercell.form_rcif_k(case_name))
    rcif_s = strip_experiment(supercell.form_rcif_supercell(case_name))
    to_super = numpy.array(case["n_cell"], dtype=float)[:, numpy.newaxis]
    n_cell_sq = float(numpy.prod(case["n_cell"])) ** 2

    intensity_k = calc_intensity(rcif_k + form_diffrn(index_hkl))
    intensity_s = calc_intensity(rcif_s + form_diffrn(index_hkl * to_super))
    # magnetic satellites are present, the index (1/4, 0, 0) is empty
    n_h = 7
    assert intensity_k[n_h:2 * n_h].max() > 0.1
    assert intensity_k[-1] == pytest.approx(0., abs=1e-10)
    numpy.testing.assert_allclose(intensity_s, n_cell_sq * intensity_k,
                                  rtol=1e-6, atol=1e-8 * intensity_s.max())

    # two domains, the second one rotated by DOMAIN_MATRIX
    intensity_kd = calc_intensity(rcif_k + form_diffrn(
        index_hkl, DOMAIN_MATRIX, DOMAIN_FRACTION))
    intensity_s2 = calc_intensity(rcif_s + form_diffrn(
        (DOMAIN_MATRIX @ index_hkl) * to_super))
    expected = (1. - DOMAIN_FRACTION) * intensity_s + DOMAIN_FRACTION * intensity_s2
    numpy.testing.assert_allclose(expected, n_cell_sq * intensity_kd,
                                  rtol=1e-6, atol=1e-8 * expected.max())
    # the domain changes the intensities (here the matrix is not a symmetry
    # of the test structure, so the nuclear ones change too)
    assert numpy.abs(intensity_kd - intensity_k)[n_h:].max() > 1e-3 * intensity_k.max()


def test_single_crystal_k_half():
    """k = (1/2, 0, 1/2) with I centring, supercell 2a x b x 2c."""
    check_with_supercell("k_half")


def test_single_crystal_k_third():
    """k = (1/3, 0, 0): both arms ±k, cos and sin, tau = 1/10, supercell 3a."""
    check_with_supercell("k_third")
