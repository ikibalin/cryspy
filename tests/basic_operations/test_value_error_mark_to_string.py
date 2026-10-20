"""value(sigma) strings written to rcif must be rounded, not truncated.

Before the fix 118.32 with sigma 0.22 was written as 118.31(22): the
fraction 118.32 % 1 * 100 = 31.999... was cut by int(), so every
read -> to_cif cycle could lower the last digit of a refined value.
"""
import pytest

from cryspy.A_functions_base.function_1_strings import (
    string_to_value_error_mark,
    value_error_mark_to_string,
)


@pytest.mark.parametrize("value, sigma, expected", [
    (118.32, 0.22, "118.32(22)"),
    (118.33, 0.22, "118.33(22)"),
    (1.15, 0.04, "1.150(40)"),
    (12.3456, 0.0029, "12.3456(29)"),
    (0.0029, 0.0029, "0.0029(29)"),
    (-0.5, 0.03, "-0.500(30)"),
    (-0.004, 0.0012, "-0.0040(12)"),
    (1234.5, 12.3, "1234(12)"),
])
def test_value_error_mark_to_string_rounds(value, sigma, expected):
    assert value_error_mark_to_string(value, sigma, "") == expected


@pytest.mark.parametrize("value, sigma", [
    (118.32, 0.22), (1.15, 0.04), (12.3456, 0.0029), (-0.004, 0.0012)])
def test_write_read_keeps_value(value, sigma):
    """Writing and reading again must not change the value: the digits
    kept are the ones of the value itself (two significant digits of
    sigma), so the round trip is exact up to float representation."""
    string = value_error_mark_to_string(value, sigma, "")
    value_read, sigma_read, _ = string_to_value_error_mark(string)
    assert value_read == pytest.approx(value, abs=1e-12)
    assert sigma_read == pytest.approx(sigma, abs=1e-12)
    # a second cycle gives the same text
    assert value_error_mark_to_string(value_read, sigma_read, "") == string
