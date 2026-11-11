"""Total moment of a position with several propagation vectors."""
import numpy
import pytest

from cryspy.D_functions_item_loop.function_1_magnetic_k import (
    calc_total_moment_range, report_total_moments)


def test_two_collinear_waves_add_and_cancel():
    """M(L) = cos(pi L_x) + cos(pi L_y) along a: 2, 0 or -2 μB."""
    amplitudes = numpy.array([[1., 1.], [0., 0.], [0., 0.]], dtype=complex)
    k_vectors = numpy.array([[0.5, 0.], [0., 0.5], [0., 0.]])
    assert calc_total_moment_range(amplitudes, k_vectors) == pytest.approx((0., 2.))


def test_two_perpendicular_waves_keep_the_size():
    """cos along a and cos along c of two k = 1/2 waves: |M| = sqrt(2) in every cell."""
    amplitudes = numpy.array([[3., 0.], [0., 0.], [0., 4.]], dtype=complex)
    k_vectors = numpy.array([[0., 0.5], [0., 0.5], [0.5, 0.]])
    assert calc_total_moment_range(amplitudes, k_vectors) == pytest.approx((5., 5.))


def test_report_groups_the_waves_by_position():
    modulation = {
        "label": ["Tb", "Tb", "Tb_2"],
        "fract_xyz": numpy.array([[0., 0., 0.5], [0.17, 0.17, 0.83], [0.25, 0.25, 0.75]]),
        "amplitude": numpy.array([[1., 0., 2.], [0., 0., 0.], [0., 1., 0.]], dtype=complex),
        "k_vector": numpy.array([[0., 0.5, 0.], [0., 0.5, 0.], [0.5, 0., 0.5]]),
    }
    text = "\n".join(report_total_moments(modulation))
    assert "| Tb | 0.0000 | 0.1700 | 0.2500 | 2 | 1.414–1.414 |" in text
    assert "| Tb_2 | 0.5000 | 0.8300 | 0.7500 | 1 | 2.000–2.000 |" in text
