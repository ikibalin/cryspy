"""Fit of the integrated intensities in the matching mode of powder patterns."""
import numpy
import pytest

from cryspy.procedure_rhochi.rhochi_pd_matching import (
    fit_matching_intensities, group_equal_positions, solve_nonnegative,
    calc_matching_signal, take_matching_mode)


def form_phase(ttheta_hkl, multiplicity, ttheta, width=0.3):
    """Gaussian peaks of unit area."""
    column = numpy.exp(-0.5 * numpy.square((ttheta[:, None] - ttheta_hkl[None, :]) / width))
    column = column / (width * numpy.sqrt(2. * numpy.pi))
    return {"column": column, "ttheta_hkl": numpy.asarray(ttheta_hkl, dtype=float),
            "multiplicity": numpy.asarray(multiplicity), "dict_in_out_phase": {}}


TTHETA = numpy.linspace(0., 20., 401)


def test_solve_nonnegative_against_unconstrained():
    rng = numpy.random.default_rng(1)
    a = rng.random((50, 4))
    x_true = numpy.array([1., 2., 0.5, 3.])
    x = solve_nonnegative(a, a @ x_true)
    assert numpy.allclose(x, x_true)
    # a negative unconstrained solution is cut to zero
    b = a @ numpy.array([1., -2., 0.5, 3.])
    x = solve_nonnegative(a, b)
    assert numpy.all(x >= 0.)
    assert x[1] == 0.
    # a column without data
    a[:, 2] = 0.
    assert solve_nonnegative(a, a @ x_true)[2] == 0.


def test_unpolarized_intensities():
    phase = form_phase(numpy.array([5., 10., 15.]), [2, 4, 6], TTHETA)
    intensity_true = numpy.array([100., 0., 30.])
    signal = phase["column"] @ intensity_true
    sigma = numpy.ones_like(signal)
    in_points = numpy.ones(signal.shape, dtype=bool)
    fit_matching_intensities([phase], signal, sigma, in_points)
    # up and down are equal halves
    assert numpy.allclose(phase["iint_plus"], 0.5 * intensity_true, atol=1e-8)
    assert numpy.allclose(phase["iint_minus"], phase["iint_plus"])
    signal_plus, signal_minus = calc_matching_signal(phase)
    assert numpy.allclose(signal_plus + signal_minus, signal)
    assert numpy.allclose(phase["dict_in_out_phase"]["iint_plus_with_factors"], phase["iint_plus"])


def test_negative_peak_gives_zero():
    phase = form_phase(numpy.array([5., 10.]), [1, 1], TTHETA)
    # a dip at 10 degrees: the intensity stays >= 0
    signal = phase["column"] @ numpy.array([50., -20.])
    in_points = numpy.ones(signal.shape, dtype=bool)
    fit_matching_intensities([phase], signal, numpy.ones_like(signal), in_points)
    assert phase["iint_plus"][1] == 0.
    assert phase["iint_plus"][0] == pytest.approx(25., rel=1e-3)


def test_excluded_points_and_reflection_outside():
    phase = form_phase(numpy.array([5., 10., 30.]), [1, 1, 1], TTHETA)
    signal = phase["column"] @ numpy.array([10., 20., 0.])
    # wrong points at 9-11 degrees are excluded: the intensity at 10 degrees
    # comes from the tails
    signal_wrong = numpy.where(numpy.abs(TTHETA - 10.) < 0.2, 1e3, signal)
    in_points = numpy.abs(TTHETA - 10.) >= 0.2
    fit_matching_intensities([phase], signal_wrong, numpy.ones_like(signal), in_points)
    assert numpy.allclose(phase["iint_plus"] * 2., [10., 20., 0.], atol=1e-6)


def test_equal_positions_share_by_multiplicity():
    phase_1 = form_phase(numpy.array([5., 5., 12.]), [2, 6, 4], TTHETA)
    phase_2 = form_phase(numpy.array([5.]), [3], TTHETA, width=0.5)
    column, index_group, share = group_equal_positions([phase_1, phase_2])
    # the reflections of one phase at the same position are merged, not the
    # reflections of another phase (another profile)
    assert column.shape[1] == 3
    assert index_group[0] == index_group[1] != index_group[3]
    assert numpy.allclose(share, [0.25, 0.75, 1., 1.])

    signal = phase_1["column"] @ numpy.array([40., 120., 10.]) + phase_2["column"] @ numpy.array([30.])
    in_points = numpy.ones(signal.shape, dtype=bool)
    fit_matching_intensities([phase_1, phase_2], signal, numpy.ones_like(signal), in_points)
    assert numpy.allclose(2. * phase_1["iint_plus"], [40., 120., 10.], atol=1e-6)
    assert numpy.allclose(2. * phase_2["iint_plus"], [30.], atol=1e-6)


def test_polarized_intensities():
    phase = form_phase(numpy.array([4., 9., 14.]), [1, 1, 1], TTHETA)
    intensity_plus = numpy.array([10., 0., 7.])
    intensity_minus = numpy.array([4., 6., 7.])
    signal_plus = phase["column"] @ intensity_plus
    signal_minus = phase["column"] @ intensity_minus
    sigma = 0.1 * numpy.ones_like(signal_plus)
    in_points = numpy.ones(signal_plus.shape, dtype=bool)
    fit_matching_intensities(
        [phase], signal_plus + signal_minus, sigma * numpy.sqrt(2.), in_points,
        signal_diff=signal_plus - signal_minus, sigma_diff=sigma * numpy.sqrt(2.))
    assert numpy.allclose(phase["iint_plus"], intensity_plus, atol=1e-6)
    assert numpy.allclose(phase["iint_minus"], intensity_minus, atol=1e-6)
    with pytest.raises(ValueError, match="sum"):
        fit_matching_intensities(
            [phase], None, None, in_points,
            signal_diff=signal_plus - signal_minus, sigma_diff=sigma)


def test_take_matching_mode():
    assert take_matching_mode({}) == "off"
    assert take_matching_mode({"matching_mode": "fixed_intensities"}) == "fixed_intensities"
    with pytest.raises(ValueError):
        take_matching_mode({"matching_mode": "le_bail"})
