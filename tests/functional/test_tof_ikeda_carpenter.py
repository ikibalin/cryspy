"""TOF powder pattern with the Ikeda-Carpenter peak shape (issue #39).

The CeCuAl model is switched to the Ikeda-Carpenter pseudo-Voigt; a
pattern calculated with known moderator parameters is fitted back from
other values, and the result survives writing to and reading from rcif.
"""
import os

import numpy
import pytest

import cryspy

DIR_TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_CECUAL = os.path.join(DIR_TESTS, "functional", "test_tof", "model_pV.rcif")

# Mantid IkedaCarpenterPV defaults (microseconds, meV)
ALPHA0, ALPHA1, BETA0, KAPPA = 1.6, 1.5, 31.9, 46.


def take_tof(obj):
    return [item for item in obj.items if isinstance(item, cryspy.TOF)][0]


def take_object():
    obj = cryspy.file_to_globaln(F_CECUAL)
    tof = take_tof(obj)
    tof.range.time_max = 5000.
    obj.fix_variables()
    profile = tof.tof_profile
    profile.peak_shape = "Ikeda-Carpenter"
    profile.alpha0, profile.alpha1 = ALPHA0, ALPHA1
    profile.beta0, profile.kappa = BETA0, KAPPA
    return obj


def calc_signal(obj):
    cryspy.rhochi_no_refinement(obj)
    tof_proc = take_tof(obj).tof_proc
    return (numpy.array(tof_proc.intensity_plus_net, dtype=float) +
            numpy.array(tof_proc.intensity_minus_net, dtype=float))


def test_ikeda_carpenter_pattern_differs_from_pseudo_voigt():
    signal_pv = calc_signal(cryspy.file_to_globaln(F_CECUAL))
    obj = take_object()
    signal_ic = calc_signal(obj)
    assert numpy.all(numpy.isfinite(signal_ic))
    assert signal_ic.max() > 0.
    n_point = signal_ic.size
    assert not numpy.allclose(signal_ic, signal_pv[:n_point])

    # rcif keeps the peak shape and its parameters
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    profile_2 = take_tof(obj_2).tof_profile
    assert profile_2.peak_shape == "Ikeda-Carpenter"
    assert float(profile_2.kappa) == pytest.approx(KAPPA)
    numpy.testing.assert_allclose(calc_signal(obj_2), signal_ic, rtol=1e-6)


def test_ikeda_carpenter_parameters_are_refined_back():
    obj = take_object()
    signal = calc_signal(obj)
    tof = take_tof(obj)
    tof_meas = tof.tof_meas
    time = numpy.array(tof.tof_proc.time, dtype=float)
    background = numpy.array(tof.tof_proc.intensity_bkg_calc, dtype=float)
    intensity = numpy.interp(
        numpy.array(tof_meas.time, dtype=float), time, signal + background)
    tof_meas.numpy_intensity = intensity
    tof_meas.numpy_intensity_sigma = numpy.sqrt(numpy.abs(intensity)) + 1.
    tof_meas.numpy_to_items()

    profile = tof.tof_profile
    profile.alpha1, profile.beta0 = 1.2 * ALPHA1, 0.8 * BETA0
    profile.alpha1_refinement, profile.beta0_refinement = True, True
    cryspy.rhochi_rietveld_refinement(obj)
    assert float(profile.alpha1) == pytest.approx(ALPHA1, rel=0.03)
    assert float(profile.beta0) == pytest.approx(BETA0, rel=0.03)
    assert float(profile.beta0_sigma) > 0.
