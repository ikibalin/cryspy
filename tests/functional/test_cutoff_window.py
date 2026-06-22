"""Focused, fast tests for the peak-range cutoff (FullProf "WDT").

They check, on small arrays, that a positive ``cutoff_fwhm``:
  * keeps the original profile shape,
  * reproduces the full profile exactly inside the window and zeros it
    outside (i.e. ``cut == full * window_mask``),
  * actually evaluates the expensive kernels on fewer points (the real
    performance win, not a post-hoc multiply).
"""
import numpy
import pytest

from cryspy.A_functions_base import powder_diffraction_tof as tof
from cryspy.A_functions_base.powder_diffraction_tof import (
    tof_Jorgensen, tof_Jorgensen_VonDreele, tof_non_convoluted_pseudo_voigt,
    calc_hpv_eta)
from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_profile_pseudo_voight, calc_h_g, calc_h_l, calc_h_pv)
from cryspy.A_functions_base.powder_diffraction_cutoff import cutoff_select

na = numpy.newaxis
SQRT_8LN2 = numpy.sqrt(8. * numpy.log(2.))


def _tof_inputs(n_t=200, n_h=12, seed=0):
    rng = numpy.random.RandomState(seed)
    time = numpy.linspace(1000., 20000., n_t)
    time_hkl = numpy.sort(rng.uniform(1500., 19500., n_h))
    alpha = rng.uniform(0.01, 0.03, n_t)
    beta = rng.uniform(0.008, 0.02, n_t)
    sigma = rng.uniform(20., 40., n_t)
    gamma = rng.uniform(10., 30., n_t)
    return time, time_hkl, alpha, beta, sigma, gamma


def _assert_equals_full_times_mask(func, mask, cutoff):
    """``func`` truncates to exactly the full result times ``mask``."""
    full = func(0.)
    cut = func(cutoff)
    assert cut.shape == full.shape
    assert mask.sum() < mask.size          # the cutoff truly removed points
    assert numpy.array_equal(cut, full * mask)


def test_tof_jorgensen_cutoff_equals_full_truncated():
    time, time_hkl, alpha, beta, sigma, _ = _tof_inputs()
    delta = time[:, na] - time_hkl[na, :]
    half_width = 3. * (sigma * SQRT_8LN2 + 1. / alpha + 1. / beta)
    mask = numpy.abs(delta) <= half_width[:, na]
    _assert_equals_full_times_mask(
        lambda cf: tof_Jorgensen(alpha, beta, sigma, time, time_hkl, cutoff_fwhm=cf),
        mask, 3.)


def test_tof_jvd_cutoff_equals_full_truncated():
    time, time_hkl, alpha, beta, sigma, gamma = _tof_inputs()
    h_pv, _ = calc_hpv_eta(sigma * SQRT_8LN2, gamma)
    delta = time[:, na] - time_hkl[na, :]
    half_width = 3. * (h_pv + 1. / alpha + 1. / beta)
    mask = numpy.abs(delta) <= half_width[:, na]
    _assert_equals_full_times_mask(
        lambda cf: tof_Jorgensen_VonDreele(
            alpha, beta, sigma, gamma, time, time_hkl, cutoff_fwhm=cf),
        mask, 3.)


def test_tof_non_conv_pv_cutoff_equals_full_truncated():
    time, time_hkl, _, _, sigma, gamma = _tof_inputs()
    h_pv, _ = calc_hpv_eta(SQRT_8LN2 * sigma, gamma)
    delta = time[:, na] - time_hkl[na, :]
    mask = numpy.abs(delta / h_pv[:, na]) <= 3.
    _assert_equals_full_times_mask(
        lambda cf: tof_non_convoluted_pseudo_voigt(sigma, gamma, time, time_hkl, cutoff_fwhm=cf),
        mask, 3.)


def test_cwl_pseudo_voigt_cutoff_equals_full_truncated():
    rng = numpy.random.RandomState(1)
    ttheta = numpy.linspace(0.1, 2.5, 200)
    ttheta_hkl = numpy.sort(rng.uniform(0.15, 2.45, 12))
    u, v, w, i_g, x, y = 0.1, -0.05, 0.05, 0.0, 0.05, 0.02
    p = (0.02, 0.01, 0.0, 0.0)
    h_g, _ = calc_h_g(u, v, w, i_g, ttheta)
    h_l, _ = calc_h_l(x, y, ttheta)
    h_pv, _ = calc_h_pv(h_g, h_l)
    delta = ttheta[:, na] - ttheta_hkl[na, :]
    z = (delta * 180. / numpy.pi) / h_pv[:, na]
    mask = numpy.abs(z) <= 30.

    def func(cf):
        return calc_profile_pseudo_voight(
            ttheta, ttheta_hkl, u, v, w, i_g, x, y, *p, cutoff_fwhm=cf)[0]

    full = func(0.)
    cut = func(30.)
    assert cut.shape == full.shape
    assert mask.sum() < mask.size
    assert numpy.array_equal(cut, full * mask)


def test_cutoff_select_gathers_in_window_subset_as_column():
    # No cutoff (callers pass inf) passes everything through unchanged.
    delta = numpy.array([[0., 5., 100.], [200., 0., 3.]])
    half_width = numpy.array([10., 10.])
    keep, params, out_delta = cutoff_select(delta, numpy.inf, half_width, (numpy.arange(2),))
    assert keep is None
    assert out_delta is delta

    # Finite cutoff: only in-window points, gathered as a column.
    keep, (gathered,), out_delta = cutoff_select(
        delta, 1., half_width, (numpy.array([10., 20.]),))
    expected = numpy.abs(delta) <= half_width[:, na]
    assert numpy.array_equal(keep, expected)
    assert out_delta.shape == (int(expected.sum()), 1)
    # per-row params gathered onto the kept points (row index of each)
    rows = numpy.nonzero(expected)[0]
    assert numpy.array_equal(gathered, numpy.array([10., 20.])[rows])


def test_tof_cutoff_kernel_work_monotonic(monkeypatch):
    # Count how many points the expensive erfc kernel is evaluated on, as
    # a timing-independent proxy for calculation cost. cutoff_fwhm = 0 must
    # evaluate the whole grid (same work as having no cutoff at all), and
    # every smaller positive cutoff must evaluate strictly fewer points
    # (so it runs faster).
    sizes = []
    orig_erfc = tof.erfc

    def spy(arg):
        sizes.append(int(numpy.size(arg)))
        return orig_erfc(arg)

    monkeypatch.setattr(tof, "erfc", spy)
    time, time_hkl, alpha, beta, sigma, gamma = _tof_inputs()
    n_full = time.size * time_hkl.size

    def kernel_work(cutoff):
        sizes.clear()
        tof_Jorgensen_VonDreele(alpha, beta, sigma, gamma, time, time_hkl, cutoff_fwhm=cutoff)
        return max(sizes)

    w_none = kernel_work(0.)   # default: no cutoff
    w_30 = kernel_work(30.)
    w_8 = kernel_work(8.)
    w_2 = kernel_work(2.)

    assert w_none == n_full     # cf=0 evaluates the full grid == no cutoff
    assert w_30 < w_none        # a cutoff reduces the work...
    assert w_8 < w_30           # ...and a smaller cutoff reduces it more
    assert w_2 < w_8
