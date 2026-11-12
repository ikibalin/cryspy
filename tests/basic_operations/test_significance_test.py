"""Hamilton's F-test of added parameters."""
import math

import pytest

from cryspy.A_functions_base.significance_test import calc_f_test_added_parameters


def test_f_value_and_p_value():
    # F = (110 - 100) / (100 / (60 - 10)) = 5; p(F(1, 50) > 5) = 0.02987 (tables)
    f_value, p_value = calc_f_test_added_parameters(110., 100., 60, 10)
    assert f_value == pytest.approx(5.)
    assert p_value == pytest.approx(0.02987, abs=2e-4)


def test_one_added_parameter_is_the_square_of_t():
    # for one parameter F = t^2 and the p-value is the two-sided one of t(n - m)
    from scipy.stats import t
    f_value, p_value = calc_f_test_added_parameters(130., 100., 25, 5)
    assert p_value == pytest.approx(2. * t.sf(math.sqrt(f_value), 20))


def test_scale_of_sigmas_does_not_matter():
    assert calc_f_test_added_parameters(110., 100., 60, 10) == pytest.approx(
        calc_f_test_added_parameters(1100., 1000., 60, 10))


def test_degenerate_cases():
    assert calc_f_test_added_parameters(100., 100., 60, 10) == (0., 1.)
    assert calc_f_test_added_parameters(100., 120., 60, 10) == (0., 1.)
    assert all(math.isnan(value) for value in calc_f_test_added_parameters(110., 100., 10, 10))
    assert calc_f_test_added_parameters(110., 0., 60, 10) == (math.inf, 0.)
