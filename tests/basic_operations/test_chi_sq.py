"""One chi-square for the refinement, the background fits and the reports."""
import numpy

from cryspy.A_functions_base.function_1_agreement_factors import (
    calc_agreement_factors, calc_chi_sq_by_signal)
from cryspy.A_functions_base.function_1_matrices import calc_chi_sq

SIGNAL_EXP = numpy.array([10., 12., numpy.nan, 30., 9.])
SIGNAL_SIGMA = numpy.array([1., 2., 1., 3., 1.5])
SIGNAL_MODEL = numpy.array([11., 11., 5., 20., 9.])
INCLUDED = numpy.array([True, True, True, True, False])


def test_value_with_excluded_and_nan_points():
    chi_sq = calc_chi_sq_by_signal(SIGNAL_EXP, SIGNAL_SIGMA, SIGNAL_MODEL,
                                   included_points=INCLUDED, flag_skip_nan=True)[0]
    assert numpy.isclose(chi_sq, 1. + 0.25 + 100./9.)
    assert numpy.isnan(calc_chi_sq_by_signal(SIGNAL_EXP, SIGNAL_SIGMA, SIGNAL_MODEL)[0])


def test_residual_max_clips_points_above_the_model():
    # residual 10/3 of the peak point counts as 1
    chi_sq = calc_chi_sq_by_signal(SIGNAL_EXP[[0, 1, 3]], SIGNAL_SIGMA[[0, 1, 3]],
                                   SIGNAL_MODEL[[0, 1, 3]], residual_max=1.)[0]
    assert numpy.isclose(chi_sq, 1. + 0.25 + 1.)


def test_derivatives_against_finite_differences():
    flag = numpy.array([True, True, False, True, True])
    signal_exp, signal_sigma = SIGNAL_EXP[flag], SIGNAL_SIGMA[flag]
    signal_model, included = SIGNAL_MODEL[flag], INCLUDED[flag]
    for residual_max in (None, 1.):
        chi_sq, dder = calc_chi_sq_by_signal(
            signal_exp, signal_sigma, signal_model, included_points=included,
            residual_max=residual_max, flag_signal_model=True)
        step = 1e-6
        for ind in range(signal_model.size):
            shift = numpy.zeros(signal_model.size)
            shift[ind] = step
            values = [calc_chi_sq_by_signal(signal_exp, signal_sigma, signal_model + sign*shift,
                                            included_points=included, residual_max=residual_max)[0]
                      for sign in (1., -1.)]
            assert numpy.isclose(dder["signal_model"][ind], (values[0] - values[1])/(2.*step),
                                 rtol=1e-6, atol=1e-6)
        assert numpy.allclose(dder["signal_model_signal_model"],
                              numpy.where(included & (dder["signal_model"] != 0.),
                                          2./numpy.square(signal_sigma), 0.))


def test_old_functions_give_the_same_chi_sq():
    flag = numpy.isfinite(SIGNAL_EXP)
    chi_sq_ref = calc_chi_sq_by_signal(SIGNAL_EXP[flag], SIGNAL_SIGMA[flag], SIGNAL_MODEL[flag])[0]
    der_y = numpy.eye(4)
    chi_sq, der_chi_sq = calc_chi_sq(SIGNAL_EXP[flag], SIGNAL_SIGMA[flag], SIGNAL_MODEL[flag], der_y)
    assert numpy.isclose(chi_sq, chi_sq_ref)
    assert numpy.allclose(der_chi_sq, -2.*(SIGNAL_EXP[flag] - SIGNAL_MODEL[flag])/numpy.square(SIGNAL_SIGMA[flag]))
    assert numpy.isclose(calc_agreement_factors(SIGNAL_EXP, SIGNAL_SIGMA, SIGNAL_MODEL)["chi_sq"], chi_sq_ref)
