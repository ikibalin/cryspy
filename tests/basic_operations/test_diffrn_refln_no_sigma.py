"""Reflections without a positive sigma are excluded from chi_sq."""
import warnings

import numpy
import pytest

from cryspy.A_functions_base.function_1_agreement_factors import calc_chi_sq_by_signal
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnReflnL

S_LOOP = """
loop_
_diffrn_refln_index_h
_diffrn_refln_index_k
_diffrn_refln_index_l
_diffrn_refln_intensity
_diffrn_refln_intensity_sigma
 0 0 2 10.0 1.0
 0 3 0  0.0 0.0
 0 0 0.5 5.0 0.5
"""


def test_reflection_with_zero_sigma_is_excluded_with_a_warning():
    loop = DiffrnReflnL.from_cif(S_LOOP)
    with pytest.warns(UserWarning, match="no positive sigma"):
        res = loop.get_dictionary()
    assert list(res["flip_ratio_excluded"]) == [False, True, False]


def test_chi_sq_is_finite_at_excluded_points_without_sigma():
    signal = numpy.array([10., 0., 5.])
    sigma = numpy.array([1., 0., 0.5])
    model = numpy.array([9., 1., 5.])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        chi_sq, dder = calc_chi_sq_by_signal(
            signal, sigma, model, included_points=numpy.array([True, False, True]),
            flag_signal_model=True)
    assert chi_sq == pytest.approx(1.)
    assert numpy.all(numpy.isfinite(dder["signal_model"]))
