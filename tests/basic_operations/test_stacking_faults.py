"""Stacking faults of layers in powder (A_functions_base/stacking_faults.py)
and the diagnosis of the peak widths (A_functions_base/peak_width_diagnostics.py)."""
import numpy
import pytest

from cryspy.A_functions_base.stacking_faults import (
    calc_fault_beta, calc_flag_faulted, calc_probability_by_width, calc_rod_geometry,
    calc_rod_mass, calc_rod_width)
from cryspy.A_functions_base.peak_width_diagnostics import (
    calc_lorentzian_columns, describe_shift, find_fault_rule, fit_sharp_and_broad,
    group_reflections)

UCP = numpy.array([10.36, 8.34, 13.39, 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])


def test_faulted_reflections_and_width():
    index_hkl = numpy.array([[0, 0, 2], [1, 1, 0], [1, 2, 1], [1, 1, 1]]).T
    shift = numpy.array([0., 0.5, 0.])
    assert list(calc_flag_faulted(index_hkl, shift)) == [False, True, False, True]
    beta = calc_fault_beta(index_hkl, shift, 0.2)
    assert beta[0] == pytest.approx(1.) and beta[1] == pytest.approx(0.6)
    d_layer = calc_rod_geometry(index_hkl, (0, 0, 2), UCP)[2]
    assert d_layer == pytest.approx(13.39 / 2.)
    gamma, delta = calc_rod_width(beta, d_layer)
    assert gamma[0] == pytest.approx(0.)
    assert gamma[1] == pytest.approx(-numpy.log(0.6) / (2. * numpy.pi * d_layer))
    assert numpy.allclose(delta, 0.)
    assert calc_probability_by_width(gamma[1], d_layer) == pytest.approx(0.2, abs=1e-4)


def test_shift_not_half_moves_the_lorentzian():
    beta = calc_fault_beta(numpy.array([[1], [0], [0]]), numpy.array([1. / 3., 0., 0.]), 0.2)
    gamma, delta = calc_rod_width(beta, 5.)
    assert gamma[0] > 0. and abs(delta[0]) > 0.


def test_rod_geometry():
    g, z, _ = calc_rod_geometry(numpy.array([[1, 1, 0], [0, 0, 2]]).T, (0, 0, 1), UCP)
    assert z[0] == pytest.approx(0.) and g[1] == pytest.approx(0., abs=1e-12)
    assert z[1] == pytest.approx(2. / 13.39)
    assert g[0] == pytest.approx(numpy.hypot(1. / 10.36, 1. / 8.34))


def test_rod_mass_keeps_the_intensity_and_matches_sampling():
    g = numpy.array([0.15, 0.15, 0.0])
    z = numpy.array([0.0, 0.1, 0.2])
    gamma = numpy.array([0.01, 0.01, 0.005])
    s_edges = numpy.linspace(0.01, 3., 30001)
    mass = calc_rod_mass(s_edges, g, z, gamma)
    assert numpy.all(mass >= -1e-12)
    assert numpy.allclose(mass.sum(axis=0), 1., atol=0.01)
    # sampling of the Lorentzian along the rod
    rng = numpy.random.default_rng(1)
    z_sample = z[1] + gamma[1] * numpy.tan(numpy.pi * (rng.random(400000) - 0.5))
    s_sample = numpy.sqrt(g[1] ** 2 + z_sample ** 2)
    hist = numpy.histogram(s_sample, bins=s_edges[::300])[0] / z_sample.size
    coarse = mass[:, 1][:30000].reshape(100, 300).sum(axis=1)
    assert numpy.abs(hist - coarse).max() < 0.005
    # a rod in the plane: sharp front at s = g, nothing below
    assert mass[s_edges[1:] < g[0] - 1e-9, 0].sum() == pytest.approx(0., abs=1e-12)


def test_narrow_lorentzian_is_a_bragg_peak():
    s_edges = numpy.linspace(0.1, 0.5, 401)
    mass = calc_rod_mass(s_edges, numpy.array([0.2]), numpy.array([0.16]), numpy.array([1e-9]))
    i_bin = numpy.searchsorted(s_edges, numpy.hypot(0.2, 0.16)) - 1
    assert mass[i_bin, 0] == pytest.approx(1., abs=1e-5)


def test_describe_shift():
    assert describe_shift(numpy.array([0., 0.5, 0.])) == "k odd"
    assert describe_shift(numpy.array([0.5, 0.5, 0.])) == "h + k odd"
    assert describe_shift(numpy.array([1. / 3., 2. / 3., 0.])) == "h - k not 3n"


def test_fit_sharp_and_broad_and_the_rule():
    x = numpy.linspace(10., 60., 2001)
    index_hkl = numpy.array([[0, 0, 2], [1, 1, 0], [1, 1, 1], [2, 0, 0], [1, 1, 2], [0, 2, 1],
                             [1, 2, 1], [1, 3, 0], [2, 0, 2]]).T
    position = numpy.linspace(15., 55., index_hkl.shape[1])
    fwhm = numpy.full(position.size, 0.3)
    flag_broad = index_hkl[1] % 2 == 1
    width_true = numpy.where(flag_broad, 4. * fwhm, fwhm)
    y = 1000. * calc_lorentzian_columns(x, position, width_true).sum(axis=1) + 50.
    sigma = numpy.sqrt(y)
    y = y + numpy.random.default_rng(2).normal(0., 1., x.size) * sigma - 50.
    column_sharp = calc_lorentzian_columns(x, position, fwhm)
    column_broad = [calc_lorentzian_columns(x, position, f * fwhm) for f in (2.5, 5.)]
    i_sharp, i_broad, s_int, _ = fit_sharp_and_broad(y, sigma, column_sharp, column_broad, x)
    part = i_broad / (i_sharp + i_broad)
    assert numpy.all(part[flag_broad] > 0.8) and numpy.all(part[~flag_broad] < 0.2)
    l_rule = find_fault_rule([index_hkl[:, [i]] for i in range(position.size)], part,
                             numpy.ones(position.size))
    assert l_rule[0]["rule"] == "k odd"
    assert numpy.allclose(l_rule[0]["shift"], (0., 0.5, 0.))


def test_group_reflections():
    group = group_reflections(numpy.array([10., 10.01, 12., 20.]), numpy.full(4, 0.2))
    assert list(group) == [0, 0, 1, 2]
