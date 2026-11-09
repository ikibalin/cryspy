"""Twinning with more than one twin law (Jana-style twin matrices)."""
import copy
import os

import numpy

import cryspy
from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")

# Twin laws as 9 components along axis 0 (row-major 3x3 matrices).
LAW_A = numpy.array([0., 0., 1., 0., 1., 0., 1., 0., 0.])  # h <-> l
LAW_B = numpy.array([-1., 0., 0., 0., 1., 0., 0., 0., -1.])  # twofold axis along b


def _dicts():
    global_dict = cryspy.load_file(F_NAME).get_dictionary()
    dict_crystal = global_dict["crystal_tbte3"]
    dict_diffrn = global_dict["diffrn_t1p8nuc"]
    return dict_diffrn, dict_crystal


def _calc(laws, fractions, flags=None):
    dict_diffrn, dict_crystal = _dicts()
    dict_diffrn = copy.deepcopy(dict_diffrn)
    dict_diffrn["twin_matrices"] = numpy.stack(laws, axis=1)
    dict_diffrn["twin_fraction"] = numpy.array(fractions, dtype=float)
    if flags is None:
        flags = [False for _ in fractions]
    dict_diffrn["flags_twin_fraction"] = numpy.array(flags, dtype=bool)
    dict_in_out = {}
    chi_sq, n_hkl, der_chi_sq, dder_chi_sq, names = calc_chi_sq_for_diffrn_by_dictionary(
        dict_diffrn, dict_crystal, dict_in_out=dict_in_out)
    ind_twin = [i for i, name in enumerate(names) if "twin_fraction" in name]
    return dict_in_out["intensity_calc"], chi_sq, der_chi_sq[ind_twin]


def test_two_twin_laws_give_weighted_sum_of_three_domains():
    i_0 = _calc([LAW_A], [0.])[0]
    i_a = _calc([LAW_A], [1.])[0]
    i_b = _calc([LAW_B], [1.])[0]
    f_a, f_b = 0.3, 0.2

    i_twin = _calc([LAW_A, LAW_B], [f_a, f_b])[0]

    i_expected = (1. - f_a - f_b) * i_0 + f_a * i_a + f_b * i_b
    assert i_twin.shape == i_0.shape
    assert numpy.allclose(i_twin, i_expected, rtol=1e-10)


def test_derivative_of_chi_sq_over_two_twin_fractions():
    fractions = numpy.array([0.3, 0.2])
    der_chi_sq = _calc([LAW_A, LAW_B], fractions, flags=[True, True])[2]
    assert der_chi_sq.shape == (2,)
    # Only the second fraction refined: its derivative must not change.
    der_chi_sq_b = _calc([LAW_A, LAW_B], fractions, flags=[False, True])[2]
    assert numpy.allclose(der_chi_sq_b, der_chi_sq[1:])

    step = 1e-6
    for i_f in range(2):
        f_plus, f_minus = fractions.copy(), fractions.copy()
        f_plus[i_f] += step
        f_minus[i_f] -= step
        chi_sq_plus = _calc([LAW_A, LAW_B], f_plus)[1]
        chi_sq_minus = _calc([LAW_A, LAW_B], f_minus)[1]
        der_numerical = (chi_sq_plus - chi_sq_minus) / (2. * step)
        assert numpy.isclose(der_chi_sq[i_f], der_numerical, rtol=1e-5)
