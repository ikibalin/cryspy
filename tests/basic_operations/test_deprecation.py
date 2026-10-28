import warnings

import pytest

import cryspy
from cryspy.A_functions_base.function_1_deprecation import deprecated


def test_deprecated_function_warns_and_returns_value():
    @deprecated
    def calc_twice(value):
        return 2 * value

    with pytest.warns(DeprecationWarning, match="calc_twice is deprecated"):
        assert calc_twice(3) == 6


def test_deprecated_method_keeps_result():
    cell = cryspy.Cell(length_a=2., length_b=2., length_c=2.,
                       angle_alpha=90., angle_beta=90., angle_gamma=90.)
    with pytest.warns(DeprecationWarning, match="Cell.calc_length_sq"):
        length_sq = cell.calc_length_sq(1., 0., 0.)
    assert length_sq == pytest.approx(4.)


def test_not_deprecated_method_is_silent():
    cell = cryspy.Cell(length_a=2., length_b=2., length_c=2.,
                       angle_alpha=90., angle_beta=90., angle_gamma=90.)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        cell.calc_sthovl(1, 0, 0)
