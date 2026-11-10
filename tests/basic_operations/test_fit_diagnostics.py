"""Signatures of misfit on synthetic Gaussian powder patterns and reflections
with known errors (shift, width, scale, background, a missing peak,
extinction)."""
import numpy

from cryspy.A_functions_base.fit_diagnostics import (
    diagnose_flip_ratios, diagnose_powder_profile, diagnose_single_crystal)

TTHETA = numpy.linspace(10., 120., 5501)
PEAKS = numpy.linspace(15., 115., 21)


def pattern(shift=0., width=0.15, scale=1., cell_strain=0., extra=None, background=100.):
    """Gaussian peaks of equal area on a flat background; cell_strain moves
    the peaks as tan(theta) (a wrong cell)."""
    y_net = numpy.zeros_like(TTHETA)
    for position in PEAKS:
        theta = numpy.radians(position) / 2.
        centre = position + shift + numpy.degrees(2. * cell_strain * numpy.tan(theta))
        y_net += scale * 2000. * numpy.exp(-0.5 * ((TTHETA - centre) / width) ** 2) / width
    if extra is not None:
        y_net += 3000. * numpy.exp(-0.5 * ((TTHETA - extra) / width) ** 2)
    return y_net, y_net + background


def diagnose(y_obs, **kwargs):
    net, y_calc = pattern()
    sigma = numpy.sqrt(numpy.maximum(y_obs, 1.))
    return dict(diagnose_powder_profile(TTHETA, y_obs, sigma, y_calc, net, **kwargs))


def test_perfect_model_has_no_signature():
    _, y_obs = pattern()
    assert set(diagnose(y_obs)) == {"FIT"}


def test_zero_shift():
    _, y_obs = pattern(shift=0.03)
    d_sign = diagnose(y_obs)
    assert "zero shift" in d_sign["SHIFT"] and "+0.03" in d_sign["SHIFT"]


def test_cell_shift_grows_with_tan_theta():
    _, y_obs = pattern(cell_strain=0.0005)
    assert "unit cell" in diagnose(y_obs)["SHIFT"]


def test_broader_peaks():
    _, y_obs = pattern(width=0.18)
    assert "broader" in diagnose(y_obs)["WIDTH"]


def test_scale_and_background():
    _, y_obs = pattern(scale=1.3, background=110.)
    d_sign = diagnose(y_obs)
    assert "1.3 times" in d_sign["SCALE"]
    assert "too low" in d_sign["BACKGROUND"]


def test_missing_peak():
    _, y_obs = pattern(extra=62.3)
    assert "62.3" in diagnose(y_obs)["UNEXPLAINED_PEAK"]


def test_intensity_trend_points_to_displacement_parameters():
    net, y_calc = pattern()
    # observed intensities fall off faster: exp(-2 B s^2) with B = 1 A^2
    sthovl = numpy.sin(numpy.radians(TTHETA) / 2.) / 1.5
    y_obs = (y_calc - 100.) * numpy.exp(-2. * sthovl ** 2) + 100.
    d_sign = diagnose(y_obs, sthovl=sthovl)
    assert "larger displacement parameters" in d_sign["INTENSITY_TREND"]


def reflections():
    rng = numpy.random.default_rng(1)
    index_hkl = rng.integers(-6, 7, size=(3, 300))
    sthovl = numpy.linalg.norm(index_hkl, axis=0) / 20. + 0.05
    i_calc = rng.uniform(1., 1000., 300)
    return index_hkl, sthovl, i_calc


def test_single_crystal_extinction():
    index_hkl, sthovl, i_calc = reflections()
    i_obs = i_calc / (1. + 0.002 * i_calc)  # secondary extinction
    i_sigma = 0.02 * i_obs + 1.
    d_sign = dict(diagnose_single_crystal(index_hkl, i_obs, i_sigma, i_calc, sthovl, n_parameter=10))
    assert "extinction" in d_sign["EXTINCTION"]
    assert "30 per parameter" in d_sign["DATA"]


def test_single_crystal_worst_reflection():
    index_hkl, sthovl, i_calc = reflections()
    i_obs = i_calc.copy()
    i_obs[7] = 3. * i_obs[7] + 100.
    i_sigma = 0.02 * i_obs + 1.
    d_sign = dict(diagnose_single_crystal(index_hkl, i_obs, i_sigma, i_calc, sthovl))
    hkl = " ".join(str(int(value)) for value in index_hkl[:, 7])
    assert d_sign["WORST"].startswith(f"reflections with the largest (I_obs - I_calc)/sigma: ({hkl})")
    assert "EXTINCTION" not in d_sign and "INTENSITY_TREND" not in d_sign


def test_flip_ratios_without_magnetic_signal():
    index_hkl = numpy.array([[1, 0, 0], [0, 1, 0], [1, 1, 1]]).transpose()
    d_sign = dict(diagnose_flip_ratios(index_hkl, numpy.array([1.2, 0.9, 1.1]),
                                       numpy.full(3, 0.02), numpy.ones(3), numpy.full(3, 0.2)))
    assert "NO_MAGNETIC_SIGNAL" in d_sign
